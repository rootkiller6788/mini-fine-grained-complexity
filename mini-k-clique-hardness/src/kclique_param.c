/**
 * kclique_param.c ? Parameterized complexity framework for k-Clique
 *
 * Implements FPT verification, kernelization bounds,
 * W[1]-hierarchy analysis, and parameterized reduction checks.
 *
 * Reference: Downey & Fellows, "Parameterized Complexity" (1999), Ch. 2, 7, 13
 * Reference: Flum & Grohe, "Parameterized Complexity Theory" (2006), Ch. 1-5
 * Reference: Chen, Huang, Kanj, Xia, "Tight lower bounds" (JCSS 2006)
 * Reference: Niedermeier, "Invitation to Fixed-Parameter Algorithms" (2006)
 *
 * Knowledge: L2 Core Concepts ? FPT, W-hierarchy, kernelization
 * Knowledge: L4 Fundamental Laws ? W[1]-completeness of k-Clique
 * Knowledge: L8 Advanced Topics ? Bounded search trees, kernel lower bounds
 */

#include "kclique_types.h"
#include "kclique_core.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ================================================================
 * FPT: Fixed-Parameter Tractability
 *
 * Definition: A parameterized problem L ? ?* ? N is Fixed-Parameter
 * Tractable (FPT) if there exists an algorithm that decides (x,k) ? L
 * in time f(k) ? |x|^{O(1)} for some computable function f.
 *
 * FPT is the class of all fixed-parameter tractable problems.
 * W[1] is the first level of the W-hierarchy (parameterized version of NP).
 *
 * Theorem (Downey & Fellows 1995): k-Clique is W[1]-complete.
 * Therefore: k-Clique ? FPT iff FPT = W[1] (which is believed false).
 * Under ETH: FPT ? W[1].
 * ================================================================ */

/**
 * Compute a lower bound on f(k): the parametric factor in any FPT
 * algorithm for k-Clique, under the assumption ETH holds.
 *
 * If k-Clique ? FPT, there would be an algorithm running in
 * time f(k) * n^c for some constant c. Under ETH, this implies
 * f(k) must grow at least as fast as 2^{?(k)}.
 *
 * Specifically: if n^{s_3 * k / log k} <= f(k) * n^c for all n,
 * then f(k) >= n^{s_3 * k / log k - c}. Since this must hold for
 * all n, f(k) must be at least 2^{?(k log k)}.
 *
 * @param k       parameter value
 * @param eth_s3  ETH constant
 * @return        minimum exponent for f(k): f(k) >= 2^{?(exponent)}
 */
double fpt_lower_bound_f_of_k(int32_t k, double eth_s3) {
    if (k <= 0 || eth_s3 <= 0.0) return 0.0;
    double logk = log((double)k + 1.0);
    if (logk < 0.01) logk = 0.01;
    /* f(k) >= 2^{s_3 * k / log k} approximately */
    return eth_s3 * (double)k / logk;
}

/**
 * Check whether a given runtime expression is FPT.
 * An algorithm with runtime O(g(n,k)) is FPT if g(n,k) = f(k) * n^{O(1)}.
 *
 * This function classifies the runtime by analyzing the separation
 * of k-dependent and n-dependent factors.
 *
 * @param n_component_exponent  exponent of n in the n-dependent part
 * @param k_component_growth    growth rate of k-dependent part
 *                              (1=exponential, 2=double-exp, etc.)
 * @return                      true if the parameterization is FPT
 */
bool is_fpt_runtime(double n_component_exponent, int32_t k_component_growth) {
    /* FPT requires: n-component is constant (independent of k),
       k-component can be anything computable but doesn't multiply
       into the exponent of n. */
    if (n_component_exponent > 10.0) return false; /* heuristic for poly */
    /* Constant exponent on n is FPT; growing with k is not */
    (void)k_component_growth; /* k-growth unrestricted in FPT */
    return true;
}

/* ================================================================
 * Kernelization
 *
 * A kernelization for a parameterized problem is a polynomial-time
 * preprocessing that reduces instance (x,k) to (x',k') where
 * |x'| <= g(k) and k' <= k, for some computable function g.
 *
 * The output (x',k') is called a kernel.
 *
 * Theorem (Bodlaender et al. 2009): k-Clique does NOT admit a
 * polynomial kernel unless NP ? coNP/poly.
 * ================================================================ */

/**
 * Check whether a given kernel size bound qualifies as a polynomial kernel.
 * A polynomial kernel has size poly(k).
 *
 * @param kernel_size  size of the kernel
 * @param k            parameter value
 * @return             true if kernel is polynomial in k
 */
bool is_polynomial_kernel(int32_t kernel_size, int32_t k) {
    if (k <= 0) return false;
    /* Polynomial kernel: size <= O(k^c) for some constant c.
       We use a generous bound of k^10. */
    double max_size = pow((double)k, 10.0);
    return (double)kernel_size <= max_size;
}

/**
 * Compute the theoretical kernel size lower bound for k-Clique.
 *
 * Theorem (Bodlaender et al. 2009, fortified by Dell & van Melkebeek 2010):
 *   k-Clique has no kernel of size O(k^{d-eps}) for any eps > 0,
 *   where d is a problem-specific constant.
 *
 * This implies k-Clique is unlikely to admit polynomial kernels,
 * which separates it from problems like k-Vertex Cover that
 * have O(k^2) kernels.
 *
 * @param k  parameter value
 * @return   minimum possible kernel size under complexity assumptions
 */
int64_t kclique_kernel_lower_bound(int32_t k) {
    if (k <= 0) return 0;
    /* Known result: no kernel of size O(k^{2-eps}) for k-Clique.
       The bound is super-polynomial. */
    return (int64_t)(pow(2.0, (double)k) / (double)k);
}

/**
 * Compare kernelization properties of k-Clique vs k-Vertex Cover.
 *
 * k-Vertex Cover: admits a 2k-vertex kernel (Nemhauser-Trotter 1975).
 * k-Clique: no polynomial kernel unless NP ? coNP/poly.
 *
 * This asymmetry shows why k-Clique is "harder" parameterized.
 *
 * @param k  parameter value
 */
void print_kernelization_comparison(int32_t k) {
    printf("Kernelization Comparison (k=%d):\n", k);
    printf("  k-Clique:       no poly kernel unless NP subset coNP/poly\n");
    printf("                  lower bound on kernel size: 2^{%d}\n", k);
    printf("  k-Vertex Cover:  O(k^2) kernel exists (2k vertices)\n");
    printf("  k-Dominating Set: no poly kernel (similar to Clique)\n");
    printf("  k-Path:          poly kernel exists\n");
}

/* ================================================================
 * W-Hierarchy and Weft
 *
 * The W-hierarchy classifies parameterized problems by the depth
 * of alternation in Boolean circuits needed to express them.
 *
 * W[1]: Problems reducible to WEIGHTED-kSAT with weft-1 circuits
 *       (only one layer of large AND gates).
 * W[2]: Weft-2 circuits (dominating set is W[2]-complete).
 * W[t]: Weft-t circuits.
 *
 * W[1] ? W[2] ? ... ? W[P] ? XP
 *
 * Theorem: k-Clique is W[1]-complete.
 * Proof: reduction from WEIGHTED-kSAT to k-Clique preserving weft.
 * ================================================================ */

/**
 * Compute the weft of a circuit that encodes k-Clique.
 *
 * The natural circuit for k-Clique has weft 2:
 * - An AND gate over all C(n,k) subsets (large fan-in, all at same level)
 * - Within each subset, AND over k*(k-1)/2 edges (large fan-in)
 * - OR gates check adjacency of each pair (small fan-in)
 *
 * However, by reorganizing, k-Clique can be expressed with weft 1
 * (which is why it's W[1]-complete, not W[2]-complete):
 * - one large AND over all pairs, nested within OR over subsets.
 *
 * @param n  number of vertices
 * @param k  clique size
 * @return   weft of the natural circuit encoding
 */
int32_t kclique_circuit_weft(int32_t n, int32_t k) {
    (void)n;
    if (k <= 1) return 0; /* trivial */
    /* k-Clique is W[1]-complete, meaning it can be expressed
       with weft-1 circuits. The standard encoding uses
       one big AND gate (weft 1) on top. */
    return 1;
}

/**
 * Classify a parameterized problem into the W-hierarchy based on
 * its circuit weft characterization.
 *
 * @param weft  circuit weft value
 * @return      string describing the W-level
 */
const char *w_hierarchy_level(int32_t weft) {
    switch (weft) {
        case 0: return "FPT (tractable)";
        case 1: return "W[1]-complete (e.g., k-Clique, k-Independent Set)";
        case 2: return "W[2]-complete (e.g., k-Dominating Set)";
        default: return "W[t]-hard (higher weft)";
    }
}

/* ================================================================
 * Parameterized Complexity of k-Clique Variants
 * ================================================================ */

/* param_result_t now defined in kclique_types.h */

/**
 * Determine parameterized complexity of k-Clique on a given graph class.
 *
 * @param max_degree    maximum degree (or -1 for unbounded)
 * @param is_planar     true if graph is planar
 * @param is_bipartite  true if graph is bipartite
 * @param has_bounded_treewidth  true if treewidth is bounded
 * @return              parameterized complexity classification
 */
param_result_t kclique_param_complexity(int32_t max_degree,
                                         bool is_planar,
                                         bool is_bipartite,
                                         bool has_bounded_treewidth) {
    /* Planar graphs: k-Clique is trivial (max clique <= 4 in planar).
       Actually: planar graphs are K_5-free, so k-Clique for k>=5 is NO.
       For k<=4: can check in O(n^4) = n^{O(1)} -> FPT (f(k)=O(1)). */
    if (is_planar) return PARAM_RESULT_FPT;

    /* Bipartite graphs: max clique <= 2. Trivial FPT. */
    if (is_bipartite) return PARAM_RESULT_FPT;

    /* Bounded treewidth: k-Clique is FPT (Courcelle's theorem + DP).
       Parameter: treewidth. Actually, k-Clique parameterized by k on
       bounded-treewidth graphs is FPT. */
    if (has_bounded_treewidth) return PARAM_RESULT_FPT;

    /* Bounded degree d: k-Clique is FPT via bounded search tree.
       Monien (1985): O(k^2 * d^k * n). This is FPT with parameter k. */
    if (max_degree > 0) return PARAM_RESULT_FPT;

    /* General graphs: W[1]-complete (Downey & Fellows 1995). */
    return PARAM_RESULT_W1_HARD;
}

/**
 * Get a human-readable description of the parameterized complexity.
 */
const char *param_result_string(param_result_t r) {
    switch (r) {
        case PARAM_RESULT_FPT:        return "FPT (fixed-parameter tractable)";
        case PARAM_RESULT_W1_HARD:    return "W[1]-hard (presumably not FPT)";
        case PARAM_RESULT_XP:         return "XP (n^{f(k)} time)";
        case PARAM_RESULT_PARA_NP:    return "para-NP-hard (NP-hard even for fixed k)";
        case PARAM_RESULT_UNKNOWN:    return "Unknown (open problem)";
        default:                      return "Invalid";
    }
}

/**
 * Print the parameterized complexity landscape for k-Clique.
 *
 * Demonstrates how structural graph parameters affect the
 * parameterized complexity of k-Clique, a key insight of
 * parameterized complexity theory.
 */
void print_kclique_param_landscape(void) {
    printf("Parameterized Complexity of k-Clique on Restricted Graph Classes\n");
    printf("===============================================================\n");
    printf("  General graphs:     W[1]-complete (Downey & Fellows 1995)\n");
    printf("  Bounded degree:     FPT (Monien 1985) - O(k^2 d^k n)\n");
    printf("  Planar graphs:      FPT (trivial, omega <= 4)\n");
    printf("  Bipartite graphs:   FPT (omega <= 2)\n");
    printf("  Bounded treewidth:  FPT (Courcelle + DP)\n");
    printf("  Interval graphs:    Polynomial (O(n+m) via perfect elimination)\n");
    printf("  Chordal graphs:     Polynomial (O(n+m) via perfect elimination)\n");
    printf("  Perfect graphs:     Polynomial (Grotschel-Lovasz-Schrijver)\n");
    printf("\n");
    printf("  Unit disk graphs:   W[1]-hard (like general graphs)\n");
    printf("  String graphs:      W[1]-hard\n");
    printf("\n");
    printf("Key insight: structural graph parameters (treewidth, degree,\n");
    printf("planarity) determine whether k-Clique is tractable or hard\n");
    printf("in the parameterized sense.\n");
}
