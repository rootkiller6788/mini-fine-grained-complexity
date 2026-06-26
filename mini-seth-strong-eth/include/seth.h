#ifndef SETH_H
#define SETH_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <math.h>

/* ============================================================================
 * mini-seth-strong-eth: Strong Exponential Time Hypothesis
 *
 * Core definitions and data structures for reasoning about the Strong
 * Exponential Time Hypothesis (SETH) and Exponential Time Hypothesis (ETH).
 *
 * Reference: Impagliazzo, Paturi, Zane (2001)
 *            "Which Problems Have Strongly Exponential Complexity?"
 *
 * SETH: For every ε > 0, there exists k such that k-SAT cannot be solved
 *       in time O(2^{(1-ε)n}) on instances with n variables.
 *
 * ETH: 3-SAT cannot be solved in time 2^{o(n)}.
 *
 * L1: Core definitions of SETH, ETH, exponential-time complexity
 * L2: Fine-grained reduction concepts, sparsification
 * L3: CNF formula data structures, variable/clause representations
 * ============================================================================ */

/* ---- L1: Core Type Definitions ---- */

/* A literal is a variable (1..n) with a sign bit */
typedef struct {
    int32_t var;       /* variable index, 1-based; 0 = end marker */
    bool    negated;   /* true if negated literal */
} literal_t;

/* A clause is a disjunction of literals */
typedef struct {
    literal_t *literals;  /* array of literals */
    int32_t    size;      /* number of literals in clause */
    bool       deleted;   /* for clause deletion during solving */
} clause_t;

/* A CNF formula: conjunction of clauses */
typedef struct {
    clause_t  *clauses;    /* array of clauses */
    int32_t    num_clauses;
    int32_t    num_vars;   /* number of variables */
    int32_t    max_clause_size;  /* k for k-SAT; 0 for general CNF */
    int32_t    capacity;   /* allocated clause capacity */
} cnf_formula_t;

/* A partial assignment: value for each variable (-1=unset, 0=F, 1=T) */
typedef struct {
    int8_t    *values;    /* values[0..num_vars-1]; -1 = unassigned */
    int32_t    num_vars;
    int32_t    num_assigned;
} assignment_t;

/* ---- L2: Exponential Time Complexity Parameters ---- */

/* Parameters characterizing the exponential complexity of k-SAT */
typedef struct {
    double    s_k;          /* inf{c : k-SAT in O(2^{c n})} */
    int32_t   k;            /* clause width */
    double    s_infinity;   /* s_∞ = lim_{k→∞} s_k */
} seth_parameters_t;

/* Result of a SAT solver run */
typedef struct {
    bool      satisfiable;
    assignment_t *witness;   /* satisfying assignment if satisfiable */
    double    time_seconds;
    uint64_t  num_branches;
    uint64_t  num_backtracks;
    uint64_t  num_propagations;
} sat_result_t;

/* ---- L1: ETH/SETH Statements ---- */

/* Check if a runtime exponent violates SETH for given k */
bool seth_violated(double achieved_exponent, int32_t k);

/* Compute the SETH limit s_k for a given k */
double seth_limit_s_k(int32_t k);

/* Check if ETH is violated: O(2^{o(n)}) for 3-SAT */
bool eth_violated(double achieved_exponent, int32_t n);

/* ---- CNF Formula Construction ---- */

cnf_formula_t *cnf_create(int32_t num_vars, int32_t capacity);
void cnf_destroy(cnf_formula_t *f);
clause_t cnf_add_clause(cnf_formula_t *f, int32_t *lits, int32_t n);
void cnf_add_clause_raw(cnf_formula_t *f, literal_t *literals, int32_t n);
cnf_formula_t *cnf_clone(const cnf_formula_t *f);
int32_t cnf_count_clauses(const cnf_formula_t *f);

/* ---- Assignment Operations ---- */

assignment_t *assign_create(int32_t num_vars);
void assign_destroy(assignment_t *a);
void assign_set(assignment_t *a, int32_t var, bool value);
void assign_unset(assignment_t *a, int32_t var);
bool assign_is_satisfying(const assignment_t *a, const cnf_formula_t *f);
int8_t assign_get(const assignment_t *a, int32_t var);
assignment_t *assign_clone(const assignment_t *a);

/* ---- Formula Simplification / Unit Propagation ---- */

/* Unit propagation: repeatedly assign forced literals.
 * Returns true if no conflict found; false if conflict detected.
 * Modifies formula in place by marking satisfied/conflicting clauses. */
bool unit_propagate(cnf_formula_t *f, assignment_t *a);

/* Pure literal elimination */
void pure_literal_eliminate(cnf_formula_t *f, assignment_t *a);

/* Check if formula is trivially satisfiable (all clauses satisfied) */
bool cnf_is_satisfied(const cnf_formula_t *f, const assignment_t *a);

/* Check if formula has a conflicting clause (all literals false) */
bool cnf_has_conflict(const cnf_formula_t *f, const assignment_t *a);

/* ---- L5: Exponential-Time SAT Algorithms ---- */

/* Brute force: O(2^n) enumeration of all assignments */
sat_result_t sat_brute_force(const cnf_formula_t *f);

/* Schöning's random walk algorithm (1999):
 * O((2(k-1)/k)^n) expected time for k-SAT
 * Reference: Schöning, "A Probabilistic Algorithm for k-SAT..."
 *            FOCS 1999 */
sat_result_t sat_schoening(const cnf_formula_t *f, int32_t max_tries,
                           int32_t max_flips);

/* PPSZ algorithm (Paturi-Pudlák-Saks-Zane 2005):
 * O(2^{(1-μ_k)n}) for k-SAT, where μ_k > 0
 * Currently best known k-SAT algorithm */
sat_result_t sat_ppsz(const cnf_formula_t *f);

/* DPLL (Davis-Putnam-Logemann-Loveland 1962):
 * Classic backtracking search with unit propagation */
sat_result_t sat_dpll(const cnf_formula_t *f);

/* CDCL (Conflict-Driven Clause Learning):
 * Modern SAT solver backbone, O(2^{n/2}) worst case
 * Reference: Marques-Silva & Sakallah (1996) */
sat_result_t sat_cdcl(const cnf_formula_t *f);

/* ---- Sparsification Lemma ---- */

/* Impagliazzo-Paturi-Zane (2001):
 * Any k-CNF formula can be reduced to O(n) linear number of
 * k-CNF formulas, each with O(n) clauses.
 * This is the key for ETH→SETH relationship. */
typedef struct {
    cnf_formula_t **formulas;
    int32_t         num_formulas;
} sparsification_result_t;

sparsification_result_t *sparsify(const cnf_formula_t *f, double epsilon);
void sparsification_free(sparsification_result_t *r);

/* ---- Isolation Lemma ---- */

/* Valiant-Vazirani (1986):
 * Any satisfiable CNF can be randomly reduced to a formula
 * with a unique satisfying assignment with prob >= 1/(4n).
 * Key for PPSZ derandomization. */
typedef struct {
    cnf_formula_t *isolated;
    int32_t        hash_bits;
    uint64_t       seed;
} isolation_result_t;

isolation_result_t *isolation_apply(const cnf_formula_t *f, uint64_t seed);
void isolation_free(isolation_result_t *ir);

/* ---- Complexity Analysis Utilities ---- */

/* Compute the running time exponent for a given algorithm */
double compute_sat_exponent(uint64_t operations, int32_t n);

/* Predict k-SAT solve time under SETH assumption */
double seth_predict_time(int32_t n, int32_t k, double base);

/* Check if a proposed algorithm would refute SETH */
bool would_refute_seth(double exponent, int32_t k, double epsilon);

/* ---- L7: Applications: Formula Statistics ---- */

double cnf_clause_density(const cnf_formula_t *f);
double ksat_phase_transition(int32_t k);
bool is_near_phase_transition(const cnf_formula_t *f);
double estimate_satisfaction_probability(const cnf_formula_t *f, int32_t num_samples);
double expected_num_solutions(int32_t n, int32_t m, int32_t k);
int32_t cnf_primal_treewidth_upper_bound(const cnf_formula_t *f);
uint64_t count_satisfying_assignments(const cnf_formula_t *f);
bool cnf_is_tautology(const cnf_formula_t *f);
double compute_sat_exponent(uint64_t operations, int32_t n);

/* ---- L8: Advanced: SETH implications ---- */

bool is_subexponential(double exponent, int32_t n);
bool is_exponential_time(double exponent, int32_t n);
bool is_linear_exponential(double exponent, int32_t n);
double seth_monotone_circuit_bound(int32_t n, int32_t k);
bool seth_implies_nc_sat_separation(int32_t k, int32_t depth);
double seth_set_partition_bound(int32_t n);

/* ---- Empirical validation ---- */

void print_seth_status_report(void);
void check_seth_consistency(
    sat_result_t (*solver)(const cnf_formula_t *),
    int32_t n, int32_t k_min, int32_t k_max, int32_t num_instances);
void print_treewidth_seth_bounds(void);
void print_complexity_dashboard(void);

#endif /* SETH_H */
