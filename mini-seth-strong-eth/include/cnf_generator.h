#ifndef CNF_GENERATOR_H
#define CNF_GENERATOR_H

#include "seth.h"

/* ============================================================================
 * CNF Formula Generators for Benchmarking and Testing
 *
 * L5: Algorithms for generating hard k-SAT instances
 * L6: Canonical problem instances for k-SAT
 * ============================================================================ */

/* ---- Random k-SAT Generators ---- */

/* Generate a random k-SAT instance with given clause-to-variable ratio.
 * Phase transition at α ≈ 4.27 for 3-SAT, α ≈ 9.93 for 4-SAT.
 * For α above threshold, formulas are almost surely unsatisfiable.
 * For α below threshold, formulas are almost surely satisfiable. */
cnf_formula_t *generate_random_ksat(int32_t num_vars, int32_t k,
                                     double clause_ratio);

/* Generate a random k-SAT instance with exactly m clauses */
cnf_formula_t *generate_random_ksat_m(int32_t num_vars, int32_t k,
                                       int32_t num_clauses);

/* ---- Structured k-SAT Instances ---- */

/* Pigeonhole principle: n+1 pigeons, n holes → unsatisfiable
 * Classic hard instance for resolution-based solvers.
 * Encoded as CNF with O(n^3) clauses. */
cnf_formula_t *generate_pigeonhole(int32_t n);

/* Parity function: XOR of all variables is true/false.
 * Encoded as CNF: requires 2^{n-1} clauses for exact representation,
 * but we encode as circuit-to-CNF for testing. */
cnf_formula_t *generate_parity_cnf(int32_t n, bool target);

/* Counting-based: exactly k out of n variables are true.
 * Cardinality constraint encoded as CNF using sequential counter. */
cnf_formula_t *generate_exactly_k(int32_t n, int32_t k);

/* Graph coloring encoded as k-colorability SAT */
cnf_formula_t *generate_graph_coloring_sat(int32_t num_vertices,
                                            int32_t num_colors,
                                            int32_t *edges, int32_t num_edges);

/* ---- Hardness-Amplifying Generators ---- */

/* Generate instances at the satisfiability threshold for given k */
cnf_formula_t *generate_threshold_ksat(int32_t num_vars, int32_t k);

/* Generate "disguised" easy instances: satisfiable but deceptive */
cnf_formula_t *generate_deceptive_ksat(int32_t num_vars, int32_t k);

/* Generate the hardest known 3-SAT instances (small n, many clauses) */
cnf_formula_t *generate_hard_3sat(int32_t num_vars);

/* ---- Utility Generators ---- */

/* Negate a CNF formula (produce DNF as CNF via fresh variables) */
cnf_formula_t *cnf_negate(const cnf_formula_t *f);

/* Expand a formula by one variable (add variable that doesn't appear) */
cnf_formula_t *cnf_expand(const cnf_formula_t *f, int32_t extra_vars);

#endif /* CNF_GENERATOR_H */
