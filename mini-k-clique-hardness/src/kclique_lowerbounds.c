/**
 * kclique_lowerbounds.c ? Detailed lower bound proofs and analysis
 *
 * Implements rigorous derivations of k-Clique lower bounds under
 * ETH/SETH, including the Chen et al. framework and Monien's
 * bounded-degree graph results.
 *
 * Reference: Chen, Huang, Kanj, Xia, "Strong computational lower
 *             bounds via parameterized complexity" (JCSS 2006)
 * Reference: Monien, "How to find long paths efficiently" (1985)
 * Reference: Lokshtanov, Marx, Saurabh, "Known algorithms on graphs
 *             of bounded treewidth are probably optimal" (SODA 2011)
 * Reference: Abboud, Backurs, Vassilevska Williams, "Tight hardness
 *             results for LCS and other sequence similarity measures"
 *             (FOCS 2015)
 *
 * Knowledge: L4 Fundamental Laws ? Detailed ETH/SETH lower bound proofs
 * Knowledge: L8 Advanced Topics ? Fine-grained reductions
 */

#include "kclique_types.h"
#include "kclique_core.h"
#include "kclique_eth.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ================================================================
 * Monien's Theorem: k-Clique in Bounded-Degree Graphs
 *
 * Theorem (Monien 1985):
 *   k-Clique in graphs with maximum degree d can be solved in
 *   O(k^2 * d^k * n) time. This is FPT with parameter k.
 *
 * Key technique: bounded search tree.
 *   - Pick a vertex v of minimum degree
 *   - Branch: either v is in the clique or not
 *   - If v is in the clique, recursively search in N(v) for (k-1)-clique
 *   - If v is NOT in the clique, remove v and search in G\{v} for k-clique
 *   - Recurrence: T(n,k) <= T(d, k-1) + T(n-1, k) + O(d + k^2)
 *   - Solving: T(n,k) <= O(k^2 * d^k * n) (FPT!)
 *
 * When d is bounded by a constant, this gives polynomial time
 * for any fixed k: O(d^k * n). But as d grows, the exponential
 * dependence on k in base d remains.
 * ================================================================ */

/**
 * Compute the runtime bound for Monien's bounded-degree k-Clique.
 *
 * T(n, k, d) = O(k^2 * d^k * n) for undirected simple graphs
 * with maximum degree d.
 *
 * @param n  number of vertices
 * @param k  clique size
 * @param d  maximum degree bound
 * @return   estimated operation count (upper bound)
 */
double monien_runtime_bound(int32_t n, int32_t k, int32_t d) {
    if (n <= 0 || k <= 0 || d <= 0) return 0.0;
    /* T(n,k,d) <= c * k^2 * d^k * n for some constant c.
       We use c=1 for the estimate. */
    double bound = (double)(k * k) * pow((double)d, (double)k) * (double)n;
    return bound;
}

/**
 * Compute the maximum degree d for which Monien's algorithm
 * outperforms the brute-force O(n^k) for given n and k.
 *
 * We want: k^2 * d^k * n < C(n,k) * k^2.
 * Simplifying: d^k * n < n^k / k!
 * => d < n / (k!)^{1/k}.
 *
 * @param n  graph size
 * @param k  clique size
 * @return   threshold degree
 */
int32_t monien_degree_threshold(int32_t n, int32_t k) {
    if (n <= 0 || k <= 0) return 0;
    if (k >= n) return 0;

    /* Compute (k!)^{1/k} approximately */
    double log_fact = 0.0;
    for (int32_t i = 1; i <= k; i++) log_fact += log((double)i);
    double kth_root = exp(log_fact / (double)k);

    double threshold = (double)n / kth_root;
    return (int32_t)threshold;
}

/* ================================================================
 * Chen et al. Lower Bound Framework
 *
 * Chen, Huang, Kanj, Xia (2006) established the canonical lower
 * bounds for k-Clique under ETH and SETH.
 *
 * Key lemma: There is an FPT reduction from kSAT to k'-Clique
 * such that:
 *   1. N = 2^{O(n/k')} vertices in the constructed graph
 *   2. k' = f(k) for some function f
 *   3. The reduction runs in time 2^{O(n/k')}
 *
 * This gives: if k'-Clique can be solved in N^{o(k')} time,
 * then kSAT can be solved in 2^{o(n)} time, contradicting ETH.
 * ================================================================ */

/**
 * Compute the instance size blowup in the SAT-to-Clique reduction.
 *
 * For a kSAT instance with n variables:
 * - Each variable appears in many clauses
 * - The reduction partitions clauses into groups
 * - Resulting graph has N = 2^{O(n/g)} vertices for some group size g
 *
 * The parameter blowup determines how finely we partition:
 * fewer groups -> larger graph but more reliable reduction.
 *
 * @param n       SAT variables
 * @param groups  number of partition groups
 * @return        number of vertices in constructed graph
 */
double chen_instance_blowup(int32_t n, int32_t groups) {
    if (n <= 0 || groups <= 0) return 0.0;
    /* N = 2^{c * n / groups} for some constant c.
       Using c=1 for estimation. */
    return exp(((double)n / (double)groups) * log(2.0));
}

/**
 * Derive the lower bound following the Chen et al. framework.
 *
 * Given:
 * - ETH: s_3 > 0
 * - SAT-to-Clique reduction with N = 2^{n/Delta} vertices
 * - k-Clique takes T(N, k) time
 *
 * If T(N, k) = N^{o(k)} = 2^{o(n * k / Delta)}, then choosing
 * Delta proportional to k gives T = 2^{o(n)}, contradicting ETH.
 *
 * Therefore: k-Clique requires N^{Omega(k)} time under ETH.
 *
 * @param n        original SAT variable count
 * @param k        clique parameter
 * @param s_3      ETH constant
 * @param delta_s  SAT-to-Clique reduction parameter
 * @return         lower bound exponent alpha
 */
double chen_eth_lower_bound(int32_t n, int32_t k,
                             double s_3, double delta_s) {
    if (n <= 0 || k <= 0 || s_3 <= 0.0 || delta_s <= 0.0) return 0.0;

    /* Under ETH, kSAT requires 2^{s_k * n} time.
       Through the reduction: 2^{s_k * n} <= (2^{n/delta_s})^{alpha * k}
       => s_k * n <= alpha * k * n / delta_s
       => alpha >= s_k * delta_s / k */

    double s_k = s_3; /* use s_3 as lower bound for s_k */
    double alpha = s_k * delta_s / (double)k;
    return alpha;
}

/* ================================================================
 * Fine-Grained Equivalence: k-Clique and Triangle Detection
 *
 * Nesetril & Poljak (1985): k-Clique is equivalent to triangle
 * detection in a suitable product graph.
 *
 * For k = 3t: k-Clique in G reduces to triangle detection in
 * G^t (the t-th tensor power of G), which has n^t vertices.
 *
 * Therefore: if triangle detection requires n^{omega-o(1)} time,
 * then k-Clique requires n^{t * omega - o(1)} = n^{omega * k/3 - o(1)}.
 *
 * This connects the fine-grained complexity of k-Clique to the
 * matrix multiplication exponent omega.
 * ================================================================ */

/**
 * Compute the Nesetril-Poljak exponent for k-Clique via
 * reduction to triangle detection.
 *
 * @param k      clique size
 * @param omega  matrix multiplication exponent
 * @return       exponent alpha such that k-Clique is in n^{alpha}
 */
double nesetril_poljak_exponent(int32_t k, double omega) {
    if (k < 3 || omega < 2.0 || omega > 3.0) return (double)k;
    /* k = 3t + r  where r in {0,1,2}
       alpha = omega * ceil(k/3) + r */
    int32_t t = (k + 2) / 3; /* ceil(k/3) */
    int32_t r = k - 3 * t + 3; /* k mod 3 with 3 mapped to 0 */
    if (r >= 3) { t++; r = 0; }
    if (r < 0) r = 0;

    return omega * (double)t + (double)r;
}

/**
 * Compare the brute-force, matrix multiplication, and SETH lower bounds
 * for small clique sizes k = 3, 4, 5, 6.
 */
void print_clique_complexity_comparison(void) {
    double omega = 2.371339; /* current best known */

    printf("k-Clique Complexity by Method\n");
    printf("====================================\n");
    printf("  k  | Brute-Force | Matrix Mult  | SETH Lower\n");
    printf("-----+-------------+--------------+-----------\n");
    for (int32_t k = 3; k <= 10; k++) {
        double brute = (double)k;
        double matrix = nesetril_poljak_exponent(k, omega);
        double seth = (double)k; /* SETH says k-o(1), essentially k */

        printf("  %2d | n^{%.1f}      | n^{%.3f}     | n^{%.1f}\n",
               k, brute, matrix, seth);
    }
    printf("\nNote: SETH says brute-force is essentially optimal.\n");
    printf("Matrix multiplication beats brute-force for small k.\n");
    printf("omega = %.6f (best known bound).\n", omega);
}

/* ================================================================
 * Conditional Lower Bound Catalog
 *
 * A comprehensive catalog of known conditional lower bounds
 * for k-Clique and related problems.
 * ================================================================ */

/**
 * Print the catalog of conditional lower bounds for k-Clique.
 *
 * This catalogs the main results in fine-grained complexity
 * theory regarding k-Clique hardness.
 */
void print_clique_lower_bound_catalog(void) {
    printf("Conditional Lower Bound Catalog for k-Clique\n");
    printf("============================================\n\n");

    printf("1. ETH => k-Clique not in n^{o(k)}\n");
    printf("   Source: Chen, Huang, Kanj, Xia (JCSS 2006)\n");
    printf("   Statement: Unless ETH fails, k-Clique has no\n");
    printf("   algorithm running in f(k) * n^{o(k)} time.\n");
    printf("   Tightness: Open (gap between lower and upper bounds)\n\n");

    printf("2. SETH => k-Clique not in n^{k-eps} for any eps>0\n");
    printf("   Source: Lokshtanov, Marx, Saurabh (SODA 2011)\n");
    printf("   Statement: Under SETH, for any constants k, eps>0,\n");
    printf("   k-Clique requires n^{k-eps} time.\n");
    printf("   Tightness: Tight! Trivial O(n^k) is optimal.\n\n");

    printf("3. W[1]-Hardness => k-Clique not in FPT\n");
    printf("   Source: Downey & Fellows (1995)\n");
    printf("   Statement: k-Clique is W[1]-complete.\n");
    printf("   Implication: No f(k)*n^{O(1)} algorithm unless W[1]=FPT.\n\n");

    printf("4. Bounded-Degree => FPT Algorithm\n");
    printf("   Source: Monien (1985)\n");
    printf("   Statement: k-Clique in degree-d graphs in O(k^2*d^k*n).\n");
    printf("   Contrast: Shows how structural parameters matter.\n\n");

    printf("5. Gap-ETH => Hardness of Approximating Clique\n");
    printf("   Source: Manurangsi-Raghavendra (FOCS 2017)\n");
    printf("   Statement: No n^{1-eps}-approximation for Max Clique\n");
    printf("   in 2^{o(n)} time under Gap-ETH.\n\n");

    printf("6. Triangle Detection Hardness\n");
    printf("   Source: Vassilevska Williams & Williams (2010)\n");
    printf("   Statement: Triangle detection and Boolean Matrix\n");
    printf("   Multiplication are subcubic-equivalent.\n\n");
}
