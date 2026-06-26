/**
 * kclique_eth.h ? ETH/SETH hardness hypothesis declarations
 *
 * Formalizes the Exponential Time Hypothesis and Strong ETH,
 * and their consequences for k-Clique lower bounds.
 *
 * Reference: Impagliazzo & Paturi, "Complexity of k-SAT" (JCSS 2001)
 * Reference: Calabro, Impagliazzo, Paturi, "Complexity of k-SAT" (2006)
 * Reference: Chen, Huang, Kanj, Xia, "Tight Lower Bounds for k-Clique" (2004)
 * Reference: Lokshtanov, Marx, Saurabh, "Lower Bounds for k-Clique" (2011)
 *
 * Knowledge: L1 Definitions ? Formal ETH and SETH statements
 * Knowledge: L2 Core Concepts ? Conditional lower bounds rationale
 * Knowledge: L4 Fundamental Laws ? ETH=>n^{Omega(k)}, SETH=>n^{k-o(1)}
 */

#ifndef KCLIQUE_ETH_H
#define KCLIQUE_ETH_H

#include <stdint.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ================================================================
 * ETH: Exponential Time Hypothesis
 *
 * Formal statement (Impagliazzo-Paturi 1999):
 *   There exists a constant s_3 > 0 such that 3SAT cannot be solved
 *   in time O(2^{s n}) for any s < s_3, where n is the number of
 *   variables. Equivalently: 3SAT requires 2^{Omega(n)} time.
 *
 * Let s_k = inf{delta : kSAT can be solved in time 2^{delta n}}.
 * ETH asserts s_3 > 0. SETH asserts lim_{k->inf} s_k = 1.
 * ================================================================ */

/**
 * Compute the ETH lower bound exponent for k-Clique.
 * Theorem (Chen et al. 2006): Under ETH, k-Clique requires
 * n^{Omega(k)} time. More precisely, for any function f(k) = o(k),
 * there is no algorithm solving k-Clique in time n^{f(k)}.
 *
 * This function computes the implied exponent given parameters.
 *
 * @param n      number of vertices in the graph
 * @param k      clique size parameter
 * @param s_3    ETH constant (lower bound for 3SAT exponent)
 * @param[out] result  pointer to store the computed lower bound value
 * @return       ETH constant s_k for the derived kSAT reduction
 */
double eth_s_k_for_clause_width(int k);

/**
 * Derive the n^{Omega(k)} lower bound from ETH.
 * Uses the Sparsification Lemma (Impagliazzo-Paturi-Zane 2001):
 * kSAT instance with n variables and m clauses can be reduced to
 * 2^{eps n} instances of kSAT with at most c(k,eps)*n clauses.
 * Then via reduction to k-Clique, we obtain the lower bound.
 *
 * Theorem (Chen et al. 2006): Unless ETH fails, k-Clique
 * cannot be solved in time f(k) * n^{o(k)} for any function f.
 *
 * @param n       graph size
 * @param k       clique size
 * @param eth_s3  assumed ETH constant s_3 (typically >= 1.0)
 * @return        the derived lower bound exponent alpha such that
 *                n^{alpha} is a lower bound up to k-dependent factors
 */
double eth_kclique_lower_bound(int32_t n, int32_t k, double eth_s3);

/**
 * Compute the exact exponent implied by ETH for a specific k.
 * For k-Clique with parameter k, ETH implies no n^{delta*k} algorithm
 * exists where delta depends on s_3. This function computes delta.
 *
 * @param k      clique size
 * @param eth_s3  ETH s_3 constant
 * @return       delta such that n^{delta*k} lower bound holds
 */
double eth_delta_exponent(int32_t k, double eth_s3);

/* ================================================================
 * SETH: Strong Exponential Time Hypothesis
 *
 * Formal statement:
 *   lim_{k->inf} s_k = 1, where s_k = inf{delta : kSAT in 2^{delta n}}.
 *
 * Equivalent formulation:
 *   For every eps > 0, there exists k such that kSAT cannot be solved
 *   in time O(2^{(1-eps)n}).
 *
 * Consequence for k-Clique (tight bound):
 *   Under SETH, for any constant k >= 3, k-Clique requires n^{k-o(1)}
 *   time. This is tight because the trivial algorithm runs in O(n^k).
 * ================================================================ */

/**
 * Compute the SETH-implied lower bound for k-Clique.
 *
 * Theorem (Chen et al. 2006, Lokshtanov et al. 2011):
 * Under SETH, k-Clique requires n^{k-o(1)} time. More precisely,
 * for every eps > 0 and every k, k-Clique cannot be solved in
 * time O(n^{k-eps}).
 *
 * @param n       graph size
 * @param k       clique size
 * @param eps     epsilon parameter from SETH (> 0)
 * @return        the lower bound exponent k - delta where delta -> 0 as n grows
 */
double seth_kclique_lower_bound(int32_t n, int32_t k, double eps);

/**
 * Analyze the tightness gap: how close is the trivial O(n^k) algorithm
 * to the SETH lower bound of n^{k-o(1)}.
 *
 * The gap is o(1) in the exponent, meaning the trivial algorithm is
 * essentially optimal under SETH.
 *
 * @param n       graph size
 * @param k       clique size
 * @return        the gap in the exponent (should approach 0 as n grows)
 */
double seth_tightness_gap(int32_t n, int32_t k);

/**
 * Derive the SAT-clause-width threshold for SETH.
 * Given eps > 0, find the minimum k such that kSAT hardness
 * implies the n^{k-eps} lower bound for k-Clique.
 *
 * @param eps     desired epsilon accuracy
 * @param n       instance size
 * @return        minimum kSAT clause width needed
 */
int32_t seth_required_sat_width(double eps, int32_t n);

/* ================================================================
 * Gap-ETH: A stronger variant
 *
 * Gap-ETH (Dinur 2016, Manurangsi-Raghavendra 2017):
 * There exists eps > 0 such that no 2^{o(n)} time algorithm can
 * distinguish between satisfiable 3SAT instances and those where
 * at most (1-eps) fraction of clauses can be simultaneously satisfied.
 *
 * Gap-ETH has strong consequences for hardness of approximation.
 * ================================================================ */

/**
 * Check if Gap-ETH implies stronger lower bounds for approximate
 * k-Clique (finding a clique of size close to the maximum).
 *
 * @param approx_factor   approximation factor gamma < 1
 * @param[out] result     resultant lower bound exponent
 * @return                true if Gap-ETH provides a meaningful bound
 */
bool gap_eth_approx_clique_bound(double approx_factor, double *result);

/**
 * Compute the parameterized complexity implication of ETH/SETH.
 * If ETH holds, then FPT != W[1]. This function evaluates the
 * evidence strength of this implication.
 *
 * @param k              parameter value
 * @param eth_s3         ETH constant
 * @param[out] fpt_bound  FPT algorithm lower bound if one exists
 * @return               confidence metric for FPT != W[1] under ETH
 */
double eth_implies_fpt_neq_w1(int32_t k, double eth_s3, double *fpt_bound);

/* ================================================================
 * Conditional lower bound: general framework
 *
 * A conditional lower bound has the form:
 *   CONJECTURE => Problem P requires Omega(n^alpha) time.
 *
 * This function encodes a specific CLB for k-Clique.
 * ================================================================ */

/**
 * Encode a conditional lower bound result as a structured object.
 *
 * @param conjecture_name  e.g., "ETH", "SETH", "3SUM-Conjecture"
 * @param target_problem   e.g., "k-Clique", "Max-Clique"
 * @param exponent         the lower bound exponent alpha
 * @param is_exponential   true if bound is 2^{Omega(n)}
 * @param is_tight         true if matching upper bound exists
 * @return                 a populated conditional_lower_bound_t structure
 *
 * Note: caller must free conjecture_name and target_problem via free_clb()
 */
#include "kclique_types.h"

conditional_lower_bound_t make_conditional_lower_bound(
    const char *conjecture_name,
    const char *target_problem,
    double exponent,
    bool is_exponential,
    bool is_tight);

/**
 * Free resources associated with a conditional lower bound.
 */
void free_clb(conditional_lower_bound_t *clb);

/**
 * Print a human-readable summary of a conditional lower bound.
 */
void print_clb(const conditional_lower_bound_t *clb);

/**
 * Compute n^alpha for a given lower bound, clamping to avoid overflow.
 * Returns -1.0 if overflow would occur.
 */
double compute_lower_bound_value(int32_t n, double alpha);

/**
 * Compare two conditional lower bounds by implied hardness.
 * Returns <0, 0, >0 for weaker, equivalent, stronger.
 */
int compare_clb_hardness(const conditional_lower_bound_t *a,
                         const conditional_lower_bound_t *b);

#ifdef __cplusplus
}
#endif

#endif /* KCLIQUE_ETH_H */
