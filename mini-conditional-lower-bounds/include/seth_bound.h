/* seth_bound.h -- SETH-based Conditional Lower Bounds
 *
 * The Strong Exponential Time Hypothesis (SETH), proposed by
 * Impagliazzo and Paturi (2001), states:
 *   For all epsilon > 0, there exists k such that k-SAT cannot
 *   be solved in time O(2^{(1-epsilon)n}).
 *
 * More formally: Let s_k = inf{delta : k-SAT can be solved in time
 * O(2^{delta n})}. SETH asserts that lim_{k->inf} s_k = 1.
 *
 * SETH is the strongest form of the Exponential Time Hypothesis (ETH),
 * which only asserts that s_3 > 0, i.e., 3-SAT requires time 2^{Omega(n)}.
 *
 * SETH implies a large number of fine-grained lower bounds:
 *   - OV requires n^{2-o(1)} (Williams, 2005)
 *   - Edit Distance requires n^{2-o(1)} (Backurs-Indyk, 2015)
 *   - LCS requires n^{2-o(1)} (Abboud-Backurs-Williams, 2015)
 *   - Frechet Distance requires n^{2-o(1)} (Bringmann, 2014)
 *   - Diameter requires n^{2-o(1)} (Roditty-Williams, 2013)
 *   - Hitting Set requires n^{2-o(1)} (Abboud et al., 2016)
 *
 * References:
 *   Impagliazzo & Paturi (2001): "On the Complexity of k-SAT",
 *        Journal of Computer and System Sciences 62(2):367-375
 *   Impagliazzo, Paturi, Zane (2001): "Which Problems Have
 *        Strongly Exponential Complexity?", JCSS 63(4):512-530
 *   Calabro, Impagliazzo, Paturi (2009): "The Complexity of
 *        Satisfiability of Small Depth Circuits"
 */

#ifndef SETH_BOUND_H
#define SETH_BOUND_H

#include "condlb.h"
#include <stdint.h>

/* ============================================================================
 * L1: Definitions - SETH and ETH formal parameters
 * ============================================================================ */

/* SAT instance: a CNF formula with n variables and m clauses */
typedef struct {
    int     n_vars;          /* number of Boolean variables */
    int     n_clauses;       /* number of clauses */
    int     max_clause_len;  /* k for k-SAT */
    int**   clauses;         /* clauses[i][j] = literal (positive=var+1, negative=-(var+1)) */
    int*    clause_sizes;    /* size of each clause */
    int     capacity;
} SatInstance;

/* SETH parameters: for a given k-SAT instance, compute s_k
 * (the infimum of delta such that k-SAT in O(2^{delta n})) */
typedef struct {
    int    k;                /* clause width */
    double s_k_inf;          /* infimum delta for this k */
    double s_k_best_upper;   /* best known upper bound delta */
    double s_k_lower_bound;  /* best unconditional lower bound delta */
    int    s_k_converging;   /* 1 if s_k appears to converge to 1 as k->inf */
} SethParameters;

/* ETH parameters: for 3-SAT specifically */
typedef struct {
    double s_3;              /* delta for 3-SAT */
    double best_upper;       /* best algorithm: O(1.307^n) or similar */
    double best_lower;       /* unconditional lower bound */
    int    eth_holds;        /* 1 if we assume ETH is true */
} EthParameters;

/* ============================================================================
 * L4: Fundamental Theorems
 * ============================================================================ */

/* Theorem (Impagliazzo-Paturi-Zane, 2001): SETH implies that CNF-SAT
 * on n variables and poly(n) clauses requires time 2^{n-o(n)}.
 *
 * This function computes the SETH lower bound for a given k.
 * Returns the conditional lower bound exponent delta_k such that
 * k-SAT requires time Omega(2^{delta_k * n}) under SETH. */
double seth_lower_bound_for_k(int k);

/* Theorem (Calabro-Impagliazzo-Paturi, 2009): SETH implies that
 * Circuit-SAT for depth-d circuits requires time
 * 2^{n - n / O(log^{d-1} n)}. */
double seth_circuit_sat_bound(int depth_d, int n_vars);

/* Theorem (Cygan et al., 2016): SETH implies lower bounds for
 * the running time of parameterized algorithms for several
 * fundamental problems. Compute the SETH lower bound for a
 * parameterized problem with parameter p. */
double seth_parameterized_bound(int n, int parameter_p);

/* ============================================================================
 * L5: Algorithms - SAT Baseline Implementations
 * ============================================================================ */

/* Create and destroy a SAT instance */
SatInstance* sat_create(int n_vars, int n_clauses, int k);
void         sat_free(SatInstance* inst);

/* Add a clause to the SAT instance. Clause is array of literals.
 * Returns clause index. */
int sat_add_clause(SatInstance* inst, const int* literals, int len);

/* Naive brute-force SAT solver: try all 2^n assignments.
 * Time: O(2^n * m). Returns 1 if satisfiable, 0 if not.
 * Also writes satisfying assignment to assignment[] if provided. */
int sat_brute_force(const SatInstance* inst, int* assignment);

/* DPLL (Davis-Putnam-Logemann-Loveland) SAT solver.
 * Uses unit propagation and pure literal elimination.
 * Much faster than brute-force in practice.
 * Time: O(2^n) worst-case, often much faster.
 * Reference: Davis, Logemann, Loveland (1962) */
int sat_dpll(const SatInstance* inst, int* assignment);

/* Schoning's randomized algorithm for k-SAT.
 * Time: O((2(k-1)/k)^n). For k=3: O(1.333^n).
 * Reference: Schoning (1999). */
int sat_schoening(const SatInstance* inst, int* assignment, int max_tries);

/* PPSZ (Paturi-Pudlak-Saks-Zane) randomized algorithm.
 * Best known worst-case bounds for k-SAT.
 * For k=3: O(1.307^n). */
int sat_ppsz(const SatInstance* inst, int* assignment);

/* Check if a given assignment satisfies a SAT instance.
 * Time: O(m * k). */
int sat_check(const SatInstance* inst, const int* assignment);

/* ============================================================================
 * L5: SETH-based Reduction Construction
 * ============================================================================ */

/* The SETH-to-OV reduction (Williams, 2005):
 * Given a k-SAT instance with n variables, construct an OV instance
 * with N = 2^{n/2} vectors of dimension d = O(m).
 * If OV can be solved in O(N^{2-epsilon}) time, then k-SAT can be
 * solved in O(2^{(1-epsilon/2)n}) time, contradicting SETH. */
void seth_to_ov(const SatInstance* sat_inst,
                int** out_vec_a, int* n_a, int* dim_a,
                int** out_vec_b, int* n_b, int* dim_b);

/* SETH to Edit Distance reduction (Backurs-Indyk, 2015):
 * Converts SAT to strings where edit distance reveals satisfiability. */
void seth_to_edit_distance(const SatInstance* sat_inst,
                           char** out_s1, int* len1,
                           char** out_s2, int* len2);

/* SETH to Hitting Set (Abboud et al., 2016):
 * SAT -> Hitting Set reduction preserving subexponential hardness. */
void seth_to_hitting_set(const SatInstance* sat_inst,
                         int** out_universe, int* universe_size,
                         int*** out_sets, int* n_sets, int** set_sizes);

/* ============================================================================
 * L7: Applications - SETH in Practice
 * ============================================================================ */

/* Verify that if an algorithm for problem P achieves time
 * O(n^{2-epsilon}), this would refute SETH.
 * Returns 1 if SETH would be refuted, 0 if consistent. */
int seth_check_refutation(double algorithm_exponent,
                          const char* problem_name,
                          const char* reduction_ref);

/* Compute the "Paturi-Zane" bound: the exact SETH threshold
 * for a given problem, i.e., the largest delta such that
 * time O(n^{delta}) is consistent with SETH. */
double seth_threshold(const char* problem_name);

/* Print a comprehensive SETH status report: current known bounds
 * for all k, convergence status, and implications. */
void seth_status_report(void);

/* ============================================================================
 * L8: Advanced Topics - Nondeterministic SETH
 * ============================================================================ */

/* NSETH (Nondeterministic SETH): Nondeterministic algorithms for
 * k-TAUT cannot achieve co-nondeterministic time 2^{(1-epsilon)n}.
 * This is the nondeterministic analog of SETH.
 * Reference: Carmosino et al. (2016). */

/* Check if a problem has a nondeterministic SETH lower bound.
 * Returns exponent delta such that nondeterministic time
 * O(2^{delta n}) is optimal under NSETH. */
double nseth_lower_bound(const char* problem_name);

/* Verify that an algorithm does NOT refute NSETH.
 * NSETH is known to be consistent with SETH but their
 * interaction is subtle. */
int nseth_consistency_check(double nondet_exponent,
                            double det_exponent);

#endif /* SETH_BOUND_H */
