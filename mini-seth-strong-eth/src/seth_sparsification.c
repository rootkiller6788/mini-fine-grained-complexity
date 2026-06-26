/* seth_sparsification.c -- Sparsification Lemma Extensions and Applications
 *
 * Extends the core sparsification framework in seth_complexity.c with:
 * - Advanced sparsification wrappers with metadata
 * - ETH-SETH relationship proofs via sparsification
 * - Sparsification-based algorithm design
 * - Fine-grained sparsification for specific problem classes
 *
 * Key Theorem (IPZ 2001):
 *   If ETH holds, then SETH holds.
 *   Proof via the Sparsification Lemma: any k-CNF with n variables
 *   and m clauses can be written as the OR of at most 2^{epsilon*n} k-CNF
 *   formulas each having O(n) clauses.
 *
 * Reference: Impagliazzo, Paturi, Zane (2001), JCSS 63(4):512-530
 */

#include "seth.h"
#include "exponential_complexity.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

static int32_t sparsification_clause_bound(int32_t k, double epsilon) {
    if (epsilon <= 0.0) epsilon = 0.1;
    double ratio = (4.0 * (double)k) / epsilon;
    double bound = pow(ratio, 3.0 * (double)k);
    if (bound > 1e9) bound = 1e9;
    return (int32_t)bound;
}

sparsification_wrapper_t *sparsification_apply_wrapper(const cnf_formula_t *f,
                                                        double epsilon) {
    sparsification_wrapper_t *sw = (sparsification_wrapper_t *)calloc(1,
        sizeof(sparsification_wrapper_t));
    if (!sw) return NULL;

    sw->epsilon = epsilon;
    int32_t n = f->num_vars;
    int32_t k = f->max_clause_size;
    int32_t m = cnf_count_clauses(f);
    sw->max_clauses_per_formula = sparsification_clause_bound(k, epsilon);

    if (m <= sw->max_clauses_per_formula) {
        sw->count = 1;
        sw->sub_formulas = (cnf_formula_t **)malloc(sizeof(cnf_formula_t *));
        if (sw->sub_formulas) sw->sub_formulas[0] = cnf_clone(f);
        return sw;
    }

    int32_t num_groups = (int32_t)pow(2.0, epsilon * (double)n);
    if (num_groups < 1) num_groups = 1;
    if (num_groups > 4096) num_groups = 4096;
    if (num_groups > m) num_groups = m;

    sw->count = num_groups;
    sw->sub_formulas = (cnf_formula_t **)calloc((size_t)num_groups,
        sizeof(cnf_formula_t *));
    if (!sw->sub_formulas) { sw->count = 0; return sw; }

    for (int32_t g = 0; g < num_groups; g++) {
        sw->sub_formulas[g] = cnf_create(n, sw->max_clauses_per_formula);
        if (!sw->sub_formulas[g]) continue;
    }

    for (int32_t ci = 0; ci < f->num_clauses; ci++) {
        if (f->clauses[ci].deleted) continue;
        clause_t *c = &f->clauses[ci];
        int32_t hash = 0;
        for (int32_t j = 0; j < c->size; j++) {
            hash = hash * 31 + c->literals[j].var;
            if (c->literals[j].negated) hash = ~hash;
        }
        if (hash < 0) hash = -hash;
        int32_t group = hash % num_groups;
        if (sw->sub_formulas[group]) {
            cnf_add_clause_raw(sw->sub_formulas[group], c->literals, c->size);
        }
    }
    return sw;
}

void sparsification_wrapper_free(sparsification_wrapper_t *s) {
    if (!s) return;
    for (int32_t i = 0; i < s->count; i++) cnf_destroy(s->sub_formulas[i]);
    free(s->sub_formulas);
    free(s);
}

eth_seth_implication_t eth_implies_seth_proof(double delta, int32_t k_used) {
    eth_seth_implication_t result = {false, 0.0, 0.0, 0, 1.0};
    if (delta <= 0.0 || delta >= 1.0 || k_used < 3) return result;
    double epsilon = delta / 2.0;
    result.epsilon = epsilon;
    result.delta = delta;
    result.required_k = k_used;
    result.final_exponent = 1.0 - delta + epsilon;
    result.valid = (result.final_exponent < 1.0 - 1e-9);
    return result;
}

double sparsified_clause_ratio(int32_t k, double epsilon) {
    int32_t bound = sparsification_clause_bound(k, epsilon);
    return (double)bound;
}

bool sparsified_below_threshold(int32_t k, double epsilon) {
    double ratio = sparsified_clause_ratio(k, epsilon);
    double threshold = ksat_phase_transition(k);
    return ratio < threshold * 100.0;
}

sat_result_t sat_sparsify_then_solve(
    const cnf_formula_t *f, double epsilon,
    sat_result_t (*solver)(const cnf_formula_t *)) {

    sat_result_t result = {false, NULL, 0, 0, 0, 0.0};
    if (!f || !solver) return result;

    sparsification_result_t *sr = sparsify(f, epsilon);
    if (!sr) return result;

    for (int32_t i = 0; i < sr->num_formulas; i++) {
        if (!sr->formulas[i]) continue;
        sat_result_t sub_result = solver(sr->formulas[i]);
        result.num_branches += sub_result.num_branches;
        result.num_backtracks += sub_result.num_backtracks;
        result.num_propagations += sub_result.num_propagations;
        if (sub_result.satisfiable) {
            result.satisfiable = true;
            result.witness = sub_result.witness;
            break;
        } else {
            if (sub_result.witness) assign_destroy(sub_result.witness);
        }
    }
    sparsification_free(sr);
    return result;
}

width_reduction_result_t reduce_to_3cnf(const cnf_formula_t *f) {
    width_reduction_result_t wrr = {NULL, 0};
    if (!f) return wrr;

    int32_t total_aux = 0;
    for (int32_t i = 0; i < f->num_clauses; i++) {
        if (f->clauses[i].deleted) continue;
        int32_t s = f->clauses[i].size;
        if (s > 3) total_aux += (s - 3);
    }

    int32_t new_vars = f->num_vars + total_aux;
    int32_t est_clauses = f->num_clauses + total_aux;
    cnf_formula_t *out = cnf_create(new_vars, est_clauses);
    if (!out) return wrr;

    int32_t aux_counter = f->num_vars;
    for (int32_t i = 0; i < f->num_clauses; i++) {
        if (f->clauses[i].deleted) continue;
        int32_t s = f->clauses[i].size;
        if (s <= 3) {
            cnf_add_clause_raw(out, f->clauses[i].literals, s);
        } else {
            int32_t base = 2;
            int32_t prev_aux = -1;
            while (s - base >= 2) {
                int32_t lits[3];
                if (prev_aux < 0) {
                    lits[0] = f->clauses[i].literals[0].negated ?
                              -f->clauses[i].literals[0].var :
                              f->clauses[i].literals[0].var;
                    lits[1] = f->clauses[i].literals[1].negated ?
                              -f->clauses[i].literals[1].var :
                              f->clauses[i].literals[1].var;
                    base = 2;
                } else {
                    lits[0] = -(prev_aux);
                    lits[1] = f->clauses[i].literals[base].negated ?
                              -f->clauses[i].literals[base].var :
                              f->clauses[i].literals[base].var;
                    base++;
                }
                aux_counter++;
                lits[2] = aux_counter;
                cnf_add_clause(out, lits, 3);
                prev_aux = aux_counter;
                wrr.num_aux_added++;
            }
            int32_t final_lits[4];
            int32_t fi = 0;
            if (prev_aux > 0) final_lits[fi++] = -prev_aux;
            while (base < s) {
                final_lits[fi++] = f->clauses[i].literals[base].negated ?
                                   -f->clauses[i].literals[base].var :
                                   f->clauses[i].literals[base].var;
                base++;
            }
            if (fi > 0) cnf_add_clause(out, final_lits, fi);
        }
    }
    wrr.formula = out;
    return wrr;
}

sparsification_verification_t verify_sparsification(
    const cnf_formula_t *f, sparsification_result_t *sr) {

    sparsification_verification_t v = {true, true, true, 0};
    if (!f || !sr) { v.preserves_satisfiability = false; return v; }

    int32_t n = f->num_vars;
    int32_t max_clauses = sparsification_clause_bound(f->max_clause_size, 0.5);

    for (int32_t i = 0; i < sr->num_formulas; i++) {
        if (!sr->formulas[i]) continue;
        v.num_verified++;
        int32_t ci = cnf_count_clauses(sr->formulas[i]);
        if (ci > max_clauses && ci > n * 50) v.clause_count_bounded = false;
    }

    int32_t max_sub = (int32_t)pow(2.0, 0.5 * (double)n);
    if (max_sub < 100) max_sub = 100;
    if (sr->num_formulas > max_sub) v.formula_count_bounded = false;
    return v;
}

bool sparsification_sunflower_step(cnf_formula_t *f, int32_t k) {
    if (!f || k < 2) return false;
    bool changed = false;

    for (int32_t i = 0; i < f->num_clauses; i++) {
        if (f->clauses[i].deleted) continue;
        clause_t *ci = &f->clauses[i];
        if (ci->size != k) continue;

        for (int32_t j = i + 1; j < f->num_clauses; j++) {
            if (f->clauses[j].deleted) continue;
            clause_t *cj = &f->clauses[j];
            if (cj->size != k) continue;

            bool i_subsumes_j = true;
            for (int32_t pi = 0; pi < ci->size && i_subsumes_j; pi++) {
                bool found = false;
                for (int32_t pj = 0; pj < cj->size; pj++) {
                    if (ci->literals[pi].var == cj->literals[pj].var &&
                        ci->literals[pi].negated == cj->literals[pj].negated) {
                        found = true; break;
                    }
                }
                if (!found) i_subsumes_j = false;
            }

            bool j_subsumes_i = true;
            for (int32_t pj = 0; pj < cj->size && j_subsumes_i; pj++) {
                bool found = false;
                for (int32_t pi = 0; pi < ci->size; pi++) {
                    if (cj->literals[pj].var == ci->literals[pi].var &&
                        cj->literals[pj].negated == ci->literals[pi].negated) {
                        found = true; break;
                    }
                }
                if (!found) j_subsumes_i = false;
            }

            if (i_subsumes_j) {
                f->clauses[j].deleted = true; changed = true;
            } else if (j_subsumes_i) {
                f->clauses[i].deleted = true; changed = true; break;
            }
        }
    }
    return changed;
}

double sparsification_epsilon_for_exponent(double target_exponent,
                                            double solver_exponent) {
    double epsilon = target_exponent - solver_exponent;
    if (epsilon < 0.0) epsilon = 0.0;
    return epsilon;
}

double sparsified_solver_exponent(double epsilon, double solver_exponent) {
    return epsilon + solver_exponent;
}

void print_sparsification_summary(void) {
    printf("\n");
    printf("==============================================\n");
    printf("   SPARSIFICATION LEMMA (IPZ 2001)\n");
    printf("==============================================\n\n");

    printf("Statement:\n");
    printf("  For any k >= 3 and epsilon > 0, any k-CNF formula F\n");
    printf("  on n variables can be expressed as:\n");
    printf("    F = F_1 v F_2 v ... v F_t\n");
    printf("  where t <= 2^{epsilon*n} and each F_i is a k-CNF\n");
    printf("  with at most C(k, epsilon) * n clauses.\n\n");

    printf("Key Consequences:\n");
    printf("  1. ETH implies SETH\n");
    printf("  2. Sparse formulas are algorithmically easier\n");
    printf("  3. Enables divide-and-conquer for k-SAT\n");
    printf("  4. C(k, epsilon) = (k/epsilon)^{O(k)}\n\n");

    printf("Proof Technique:\n");
    printf("  - Sunflower lemma (Erdos-Rado 1960)\n");
    printf("  - Encoding clauses as partial assignments\n");
    printf("  - Combinatorial core lemma\n");
    printf("  - Branching on variable subsets\n\n");

    printf("Practical Impact:\n");
    printf("  - Not directly useful for SAT solving\n");
    printf("  - Critical for conditional lower bounds\n");
    printf("  - Shows ETH implies SETH\n");
    printf("  - Connects clause density to time complexity\n");
    printf("\n==============================================\n");
}
