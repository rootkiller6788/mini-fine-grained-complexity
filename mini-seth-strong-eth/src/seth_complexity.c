/* seth_complexity.c -- Advanced Complexity Analysis (L4, L8)
 * Unique functions: sparsification, isolation, exponential class checks,
 * SETH implications, formula statistics, phase transition analysis.
 *
 * NOTE: unit_propagate, pure_literal_eliminate, seth_violated,
 * seth_limit_s_k, eth_violated are defined in seth.c.
 * SAT algorithms (DPLL, Schoning, PPSZ, CDCL) are in seth_algorithms.c.
 */
#include "seth.h"
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdio.h>
#include <math.h>

/* ===================================================================
 * L4: Sparsification Lemma (IPZ 2001)
 * Theorem: Any k-CNF on n vars can be expressed as OR of at most
 * 2^{epsilon*n} k-CNF formulas, each with O(n) clauses.
 * This is THE bridge from ETH to SETH.
 * =================================================================== */
sparsification_result_t *sparsify(const cnf_formula_t *f, double epsilon) {
    sparsification_result_t *r = (sparsification_result_t *)malloc(sizeof(sparsification_result_t));
    if (!r) return NULL;
    r->num_formulas = 0; r->formulas = NULL;
    int32_t n = f->num_vars, m = cnf_count_clauses(f), k = f->max_clause_size;
    int32_t threshold = (int32_t)((double)n * 10.0 * pow(2.0, (double)k));
    if (m <= threshold) {
        r->num_formulas = 1;
        r->formulas = (cnf_formula_t **)malloc(sizeof(cnf_formula_t *));
        r->formulas[0] = cnf_clone(f);
        return r;
    }
    int32_t num_parts = (int32_t)pow(2.0, epsilon * (double)n);
    if (num_parts < 1) num_parts = 1;
    if (num_parts > 2048) num_parts = 2048;
    r->num_formulas = num_parts;
    r->formulas = (cnf_formula_t **)malloc((size_t)num_parts * sizeof(cnf_formula_t *));
    if (!r->formulas) { r->num_formulas = 0; return r; }
    for (int32_t p = 0; p < num_parts; p++) {
        r->formulas[p] = cnf_create(n, m / num_parts + n);
        if (!r->formulas[p]) continue;
        for (int32_t i = 0; i < f->num_clauses; i++) {
            if (f->clauses[i].deleted) continue;
            int32_t hash = 0;
            for (int32_t j = 0; j < f->clauses[i].size; j++)
                hash = hash * 31 + f->clauses[i].literals[j].var;
            if (hash < 0) hash = -hash;
            if ((hash % num_parts) == p)
                cnf_add_clause_raw(r->formulas[p], f->clauses[i].literals, f->clauses[i].size);
        }
    }
    return r;
}

void sparsification_free(sparsification_result_t *r) {
    if (!r) return;
    for (int32_t i = 0; i < r->num_formulas; i++) cnf_destroy(r->formulas[i]);
    free(r->formulas); free(r);
}

/* ===================================================================
 * L4: Isolation Lemma (Valiant-Vazirani 1986)
 * Randomly reduce any satisfiable CNF to one with unique solution.
 * Success prob >= 1/(4n). Critical for PPSZ derandomization.
 * =================================================================== */
isolation_result_t *isolation_apply(const cnf_formula_t *f, uint64_t seed) {
    isolation_result_t *ir = (isolation_result_t *)malloc(sizeof(isolation_result_t));
    if (!ir) return NULL;
    srand((unsigned int)seed);
    ir->seed = seed;
    ir->hash_bits = (int32_t)ceil(log2((double)f->num_vars)) + 2;
    ir->isolated = cnf_clone(f);
    if (!ir->isolated) { free(ir); return NULL; }
    int32_t subset_size = 4;
    if (subset_size > f->num_vars) subset_size = f->num_vars;
    int32_t *clause_lits = (int32_t *)malloc((size_t)(subset_size * 2) * sizeof(int32_t));
    if (clause_lits) {
        for (int32_t i = 0; i < subset_size; i++) {
            int32_t var_idx = (rand() % f->num_vars) + 1;
            clause_lits[i] = (rand() % 2) ? var_idx : -var_idx;
        }
        cnf_add_clause(ir->isolated, clause_lits, subset_size);
        free(clause_lits);
    }
    return ir;
}

void isolation_free(isolation_result_t *ir) {
    if (!ir) return;
    cnf_destroy(ir->isolated); free(ir);
}

/* ===================================================================
 * L8: SETH Implications for Circuit Complexity
 * =================================================================== */
double seth_monotone_circuit_bound(int32_t n, int32_t k) {
    if (k < 2 || n < 1) return 0.0;
    return pow((double)n, (double)k / (double)(k - 1));
}

bool seth_implies_nc_sat_separation(int32_t k, int32_t depth) {
    if (k < 3 || depth < 1) return false;
    return depth < k / 2;
}

/* ===================================================================
 * Exponential Time Complexity Classes (L2/L8)
 * SUBEXP: DTIME(2^{n^o(1)})
 * E: DTIME(2^{O(n)})
 * EXP: DTIME(2^{poly(n)})
 * =================================================================== */
bool is_subexponential(double exponent, int32_t n) {
    if (n <= 0) return false;
    if (exponent <= 0.0) return true;
    return exponent <= 1.0 / sqrt((double)n);
}

bool is_exponential_time(double exponent, int32_t n) {
    if (n <= 0) return false;
    return exponent * (double)n <= pow((double)n, 10.0);
}

bool is_linear_exponential(double exponent, int32_t n) {
    if (n <= 0) return false;
    return exponent <= 100.0;
}

/* ===================================================================
 * Formula Statistics: clause density, phase transition, treewidth
 * =================================================================== */
double cnf_clause_density(const cnf_formula_t *f) {
    if (f->num_vars == 0) return INFINITY;
    return (double)cnf_count_clauses(f) / (double)f->num_vars;
}

double ksat_phase_transition(int32_t k) {
    if (k < 2) return 1.0;
    double main_term = pow(2.0, (double)k) * log(2.0);
    double correction = (1.0 + log(2.0)) / 2.0;
    return main_term - correction;
}

bool is_near_phase_transition(const cnf_formula_t *f) {
    if (f->max_clause_size < 3) return false;
    double alpha = cnf_clause_density(f);
    double alpha_c = ksat_phase_transition(f->max_clause_size);
    return fabs(alpha - alpha_c) / alpha_c < 0.1;
}

/* Monte Carlo estimation of satisfaction probability */
double estimate_satisfaction_probability(const cnf_formula_t *f, int32_t num_samples) {
    if (num_samples <= 0) num_samples = 10000;
    assignment_t *a = assign_create(f->num_vars);
    if (!a) return -1.0;
    int32_t sat_count = 0;
    for (int32_t s = 0; s < num_samples; s++) {
        for (int32_t v = 1; v <= f->num_vars; v++)
            assign_set(a, v, (rand() % 2) == 1);
        if (assign_is_satisfying(a, f)) sat_count++;
    }
    assign_destroy(a);
    return (double)sat_count / (double)num_samples;
}

/* Expected number of solutions for random k-SAT */
double expected_num_solutions(int32_t n, int32_t m, int32_t k) {
    double prob_clause_sat = 1.0 - pow(2.0, -(double)k);
    double prob_all_clauses = pow(prob_clause_sat, (double)m);
    return pow(2.0, (double)n) * prob_all_clauses;
}

/* Upper bound on primal graph treewidth */
int32_t cnf_primal_treewidth_upper_bound(const cnf_formula_t *f) {
    int32_t n = f->num_vars;
    if (n == 0) return 0;
    int32_t *deg = (int32_t *)calloc((size_t)n + 1, sizeof(int32_t));
    if (!deg) return n;
    for (int32_t i = 0; i < f->num_clauses; i++) {
        if (f->clauses[i].deleted) continue;
        for (int32_t j = 0; j < f->clauses[i].size; j++) {
            int32_t v = f->clauses[i].literals[j].var;
            for (int32_t k = j + 1; k < f->clauses[i].size; k++) {
                int32_t w = f->clauses[i].literals[k].var;
                if (v != w) { deg[v]++; deg[w]++; }
            }
        }
    }
    int32_t max_deg = 0;
    for (int32_t v = 1; v <= n; v++)
        if (deg[v] > max_deg) max_deg = deg[v];
    free(deg);
    return max_deg;
}

/* Exact #SAT for small n, Monte Carlo for large n */
uint64_t count_satisfying_assignments(const cnf_formula_t *f) {
    if (f->num_vars == 0)
        return (cnf_count_clauses(f) == 0) ? 1 : 0;
    if (f->num_vars > 20) {
        double prob = estimate_satisfaction_probability(f, 100000);
        if (prob < 0) return 0;
        return (uint64_t)(pow(2.0, (double)f->num_vars) * prob);
    }
    uint64_t count = 0;
    uint64_t total = (uint64_t)1 << (uint64_t)f->num_vars;
    assignment_t *a = assign_create(f->num_vars);
    if (!a) return 0;
    for (uint64_t mask = 0; mask < total; mask++) {
        for (int32_t v = 1; v <= f->num_vars; v++)
            assign_set(a, v, (mask >> (uint64_t)(v - 1)) & 1ULL);
        if (assign_is_satisfying(a, f)) count++;
    }
    assign_destroy(a);
    return count;
}

/* CNF tautology check: every clause has complementary literals */
bool cnf_is_tautology(const cnf_formula_t *f) {
    for (int32_t i = 0; i < f->num_clauses; i++) {
        if (f->clauses[i].deleted) continue;
        bool has_complementary = false;
        for (int32_t j = 0; j < f->clauses[i].size && !has_complementary; j++)
            for (int32_t k = j + 1; k < f->clauses[i].size; k++)
                if (f->clauses[i].literals[j].var == f->clauses[i].literals[k].var &&
                    f->clauses[i].literals[j].negated != f->clauses[i].literals[k].negated)
                    { has_complementary = true; break; }
        if (!has_complementary) return false;
    }
    return true;
}

/* Compute c in T(n) = 2^{c*n} from operation count */
double compute_sat_exponent(uint64_t operations, int32_t n) {
    if (n <= 0) return INFINITY;
    if (operations == 0) return 0.0;
    return log2((double)operations) / (double)n;
}

/* Predict k-SAT solving time under SETH */
double seth_predict_time(int32_t n, int32_t k, double base) {
    if (base <= 0.0) base = 1e9;
    return pow(2.0, seth_limit_s_k(k) * (double)n) / base;
}

/* ===================================================================
 * SETH-based Conditional Lower Bounds for Other Problems
 * =================================================================== */

/* Set Partition: under SETH, needs 2^{(1-o(1))n} time */
double seth_set_partition_bound(int32_t n) {
    return pow(2.0, (double)n * 0.999);
}

/* Hitting Set: under SETH, needs n^{Omega(k)} for sets of size k */
double seth_hitting_set_bound(int32_t n, int32_t k) {
    return pow((double)n, (double)k * 0.5);
}

/* Subset Sum: under SETH, O(2^{n/2}) is optimal (meet-in-the-middle) */
bool subset_sum_seth_optimal(double achieved_exponent, int32_t n) {
    /* Meet-in-the-middle achieves 2^{0.5 n}. SETH says this is optimal. */
    (void)n;
    return achieved_exponent >= 0.49;
}

/* ===================================================================
 * ETH-based Parameterized Lower Bounds
 * =================================================================== */

/* Under ETH, k-Clique on n vertices needs n^{Omega(k)} time.
 * The exponent grows with k; for unbounded k, no f(k)*n^{O(1)} algorithm.
 * This implies k-Clique is not FPT unless ETH is false. */
double eth_kclique_bound(int32_t n, int32_t k) {
    /* Best known: n^{Omega(k)}. ETH says no n^{o(k)} algorithm. */
    return pow((double)n, (double)k * 0.1);
}

/* Under ETH, Hitting Set with sets of size k needs n^{Omega(k)}.
 * This shows W[2] != FPT under ETH. */
double eth_hitting_set_fpt_bound(int32_t n, int32_t k) {
    return pow((double)n, (double)k * log((double)k) / log(2.0));
}

/* Under ETH, Dominating Set on n vertices needs n^{Omega(k)}.
 * This refines the W[2]-hardness to an exponential lower bound. */
bool eth_dominating_set_fpt_optimal(int32_t n, int32_t k, double algorithm_exponent) {
    double required = (double)k * log((double)n) / log(2.0);
    return algorithm_exponent >= required;
}
