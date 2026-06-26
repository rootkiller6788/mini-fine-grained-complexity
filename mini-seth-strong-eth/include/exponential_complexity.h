#ifndef EXPONENTIAL_COMPLEXITY_H
#define EXPONENTIAL_COMPLEXITY_H

#include "seth.h"

/* ============================================================================
 * Exponential Complexity Analysis Framework
 *
 * Additional specialized types and functions for exponential-time
 * complexity analysis. Core types are in seth.h.
 *
 * L4: Fundamental laws connecting SETH/ETH to algorithm lower bounds
 * L8: Advanced topics in exponential-time complexity
 *
 * Reference: Impagliazzo, Paturi (1999) - "Complexity of k-SAT"
 *            Calabro, Impagliazzo, Paturi (2006) - "Duality between
 *            clause width and time for SAT"
 * ============================================================================ */

/* ---- Sparsification Wrapper (distinct from sparsification_result_t) ---- */

/* Higher-level sparsification with metadata tracking */
typedef struct {
    cnf_formula_t **sub_formulas;
    int32_t         count;
    double          epsilon;
    int32_t         max_clauses_per_formula;
} sparsification_wrapper_t;

sparsification_wrapper_t *sparsification_apply_wrapper(const cnf_formula_t *f,
                                                        double epsilon);
void sparsification_wrapper_free(sparsification_wrapper_t *s);

/* ---- ETH-to-SETH Implication Proof ---- */

typedef struct {
    bool     valid;
    double   epsilon;
    double   delta;
    int32_t  required_k;
    double   final_exponent;
} eth_seth_implication_t;

/* Validate the ETH-implies-SETH proof structure for given parameters */
eth_seth_implication_t eth_implies_seth_proof(double delta, int32_t k_used);

/* ---- Width Reduction: General CNF to 3-CNF ---- */

typedef struct {
    cnf_formula_t *formula;
    int32_t         num_aux_added;
} width_reduction_result_t;

/* Reduce a general CNF to equisatisfiable 3-CNF using Tseitin encoding */
width_reduction_result_t reduce_to_3cnf(const cnf_formula_t *f);

/* ---- Sparsification Verification ---- */

typedef struct {
    bool preserves_satisfiability;
    bool clause_count_bounded;
    bool formula_count_bounded;
    int32_t num_verified;
} sparsification_verification_t;

/* Verify that sparsification preserves satisfiability and bounds */
sparsification_verification_t verify_sparsification(
    const cnf_formula_t *f, sparsification_result_t *sr);

/* ---- Sunflower-Based Sparsification ---- */

/* Apply sunflower collapsing step to reduce clause count */
bool sparsification_sunflower_step(cnf_formula_t *f, int32_t k);

/* ---- Clause-Variable Ratio Analysis ---- */

double sparsified_clause_ratio(int32_t k, double epsilon);
bool sparsified_below_threshold(int32_t k, double epsilon);

/* ---- Sparsification-Driven Algorithm ---- */

sat_result_t sat_sparsify_then_solve(
    const cnf_formula_t *f, double epsilon,
    sat_result_t (*solver)(const cnf_formula_t *));

/* ---- Sparsification Prediction ---- */

double sparsification_epsilon_for_exponent(double target_exponent,
                                            double solver_exponent);
double sparsified_solver_exponent(double epsilon, double solver_exponent);

/* ---- Print Functions ---- */

/* Print a detailed summary of the sparsification lemma */
void print_sparsification_summary(void);

/* ---- ETH Parameter Report ---- */

/* Print a report of current ETH/SETH status */
void print_complexity_dashboard(void);

#endif /* EXPONENTIAL_COMPLEXITY_H */
