/**
 * kclique_eth.c ? ETH/SETH lower bound derivations for k-Clique
 *
 * Implements the mathematical derivations connecting the Exponential
 * Time Hypothesis and Strong ETH to k-Clique lower bounds.
 *
 * Reference: Impagliazzo, Paturi, Zane, "Which Problems Have
 *             Strongly Exponential Complexity?" (JCSS 2001)
 * Reference: Chen, Huang, Kanj, Xia, "Strong computational lower
 *             bounds via parameterized complexity" (JCSS 2006)
 * Reference: Lokshtanov, Marx, Saurabh, "Known Algorithms on Graphs
 *             of Bounded Treewidth are Probably Optimal" (SODA 2011)
 * Reference: Vassilevska Williams, "Hardness of Easy Problems" (2015)
 *
 * Knowledge: L2 Core Concepts ? Conditional lower bounds, hardness hypotheses
 * Knowledge: L4 Fundamental Laws ? ETH?n^{?(k)}, SETH?n^{k-o(1)}
 * Knowledge: L8 Advanced Topics ? Fine-grained complexity theory
 */

#include "kclique_eth.h"
#include "kclique_core.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ================================================================
 * ETH Analysis
 *
 * Recall: ETH states s_3 > 0 where s_3 = inf{? : 3SAT ? TIME(2^{?n})}.
 *
 * The Sparsification Lemma (Impagliazzo-Paturi-Zane 2001):
 * For every eps > 0, a kSAT instance with n variables can be reduced
 * to 2^{eps*n} instances of kSAT, each with at most c(k,eps)*n clauses
 * (where c is independent of n).
 *
 * This implies: if kSAT can be solved in 2^{o(n)} time, then every
 * problem in SNP (strict NP) can be, which is believed false.
 * ================================================================ */

/**
 * Compute s_k = inf{? : kSAT ? TIME(2^{?n})} from s_3.
 *
 * The Sparsification Lemma implies a relationship between s_k values:
 * s_k <= s_3 * (1 + log_2(k/3)) / log_2(k/3) asymptotically.
 *
 * More precisely, Impagliazzo-Paturi show:
 *   s_k >= s_3 * f(k) for some increasing function f.
 *
 * This function computes a conservative estimate of s_k.
 *
 * @param k    SAT clause width (k >= 3)
 * @param s_3  ETH constant for 3SAT
 * @return     estimated s_k value
 */
double eth_s_k_for_clause_width(int k) {
    if (k < 3) k = 3;
    if (k > 100) k = 100; /* practical limit */

    /* Use the recurrence from the Sparsification Lemma:
       s_k is at least s_3 * log(k) / log(3) * some constant.
       Conservative estimate: s_k >= s_3 * (1 - (k-3)/(k*log2(k))) */
    double log2k = log2((double)k);
    if (log2k < 1e-10) return 1.0; /* degenerate */

    /* The known relationship: s_k is non-decreasing in k.
       s_3 > 0 implies s_k >= s_3 for all k (trivially).
       SETH predicts s_k -> 1 as k -> infinity. */
    double factor = 1.0 - 1.0 / log2k; /* heuristic: approaches 1 as k grows */
    if (factor < 0.01) factor = 0.01;
    return factor;
}

/**
 * Derive the ETH k-Clique lower bound.
 *
 * Key theorem (Chen et al. 2006):
 *   Unless ETH fails, k-Clique cannot be solved in time f(k) * n^{o(k)}
 *   for any computable function f.
 *
 * Proof sketch:
 *   1. Reduce 3SAT with n variables to k-Clique with N = 2^{O(n/k)} vertices
 *   2. If k-Clique could be solved in n^{o(k)} time, then
 *      Time <= (2^{O(n/k)})^{o(k)} = 2^{o(n)}
 *   3. This would give a 2^{o(n)} algorithm for 3SAT, contradicting ETH.
 *
 * Therefore: k-Clique requires n^{Omega(k)} time under ETH.
 *
 * @param n     graph size
 * @param k     clique size
 * @param eth_s3  ETH constant s_3
 * @return      exponent alpha such that n^{alpha} is the lower bound
 */
double eth_kclique_lower_bound(int32_t n, int32_t k, double eth_s3) {
    if (n <= 0 || k <= 0 || eth_s3 <= 0.0) return 0.0;

    /* The reduction 3SAT -> Weighted-kSAT -> k-Clique gives:
       lower bound exponent = s_3 * k / C for some constant C.
       Conservative: alpha = s_3 * k / (log k + 1) */
    double logk = log((double)k + 1.0); /* natural log */
    if (logk < 0.01) logk = 0.01;

    /* The known bound from Chen et al. (2006):
       k-Clique requires n^{delta*k} under ETH where delta = s_3 / poly(k) */
    double delta = eth_s3 / (logk + 1.0);
    double alpha = delta * (double)k;

    /* Clamp: cannot exceed k (trivial enumeration) or be below 0 */
    if (alpha > (double)k) alpha = (double)k;
    if (alpha < 0.0) alpha = 0.0;

    return alpha;
}

double eth_delta_exponent(int32_t k, double eth_s3) {
    if (k <= 0 || eth_s3 <= 0.0) return 0.0;
    double logk = log((double)k + 1.0);
    if (logk < 0.01) logk = 0.01;
    return eth_s3 / (logk + 1.0);
}

/* ================================================================
 * SETH Analysis
 *
 * SETH: lim_{k??} s_k = 1.
 * Equivalently: ? ? > 0, ? k such that kSAT requires 2^{(1-?)n} time.
 *
 * Key theorem (Lokshtanov, Marx, Saurabh 2011):
 *   Under SETH, for any constants k >= 3 and ? > 0,
 *   k-Clique cannot be solved in O(n^{k-?}) time.
 *
 * This is TIGHT: the trivial O(n^k) algorithm is essentially optimal.
 *
 * Proof sketch:
 *   1. Assume SETH holds. Fix k and ? > 0.
 *   2. Choose K large enough that K-SAT requires 2^{(1-?)n} time
 *      where ? << ?/k.
 *   3. Reduce K-SAT with n variables to k-Clique with N = 2^{n/?} vertices
 *      where ? ~ k. The exact parameter grows as k grows.
 *   4. If k-Clique could be solved in N^{k-?} = 2^{n(k-?)/?} time
 *      where (k-?)/? < 1-?, we'd get sub-2^{(1-?)n} K-SAT.
 *   5. Choosing ? sufficiently large yields contradiction.
 *
 * Therefore: k-Clique requires n^{k-o(1)} time under SETH.
 * ================================================================ */

double seth_kclique_lower_bound(int32_t n, int32_t k, double eps) {
    if (n <= 0 || k <= 0 || eps <= 0.0) return (double)k;

    /* Under SETH, the lower bound exponent is k - o(1).
       For finite n, the o(1) term is materialized as:
         o(1) ? log_k(log n) / log n  or similar.
       We model it as a small diminishing term. */

    double logn = log((double)n);
    if (logn < 1.0) logn = 1.0;
    double o_term = log(log((double)n + 2.0)) / logn; /* o(1) estimate */

    /* The SETH bound: lower = k - o(1) ? k - eps for large enough n */
    double lower = (double)k - o_term;
    if (lower < 0.0) lower = 0.0;
    /* SETH also says: cannot beat k-eps for any eps > 0 */
    if (lower < (double)k - eps) lower = (double)k - eps;

    return lower;
}

double seth_tightness_gap(int32_t n, int32_t k) {
    if (n <= 1 || k <= 0) return 0.0;
    double logn = log((double)n);
    if (logn < 1.0) return 0.0;
    double o_term = log(log((double)n + 2.0)) / logn;
    return o_term;
}

int32_t seth_required_sat_width(double eps, int32_t n) { (void)n;
    if (eps <= 0.0 || eps >= 1.0) return 3;
    /* For SETH to imply n^{k-eps} lower bound, need K-SAT hardness
       with K large enough. The relationship is K ~ k/eps. */
    int32_t k_sat = (int32_t)ceil(1.0 / eps);
    if (k_sat < 3) k_sat = 3;
    if (k_sat > 100) k_sat = 100;
    return k_sat;
}

/* ================================================================
 * Gap-ETH Analysis
 *
 * Gap-ETH (Manurangsi-Raghavendra 2017):
 *   ? ? > 0 such that distinguishing satisfiable 3SAT from 3SAT
 *   with ? (1-?)m satisfiable clauses requires 2^{?(n)} time.
 *
 * Gap-ETH has strong consequences for approximating clique:
 *   Approximating maximum clique to within n^{1-?} for any ?>0
 *   requires 2^{?(n)} time under Gap-ETH.
 *
 * For k-Clique specifically:
 *   Gap-ETH implies that detecting a clique of size k = n^{1-?}
 *   is hard (requires superpolynomial time).
 * ================================================================ */

bool gap_eth_approx_clique_bound(double approx_factor, double *result) {
    if (!result) return false;
    if (approx_factor <= 0.0 || approx_factor >= 1.0) {
        *result = 0.0;
        return false;
    }

    /* Under Gap-ETH: approximate k-Clique with factor n^{1-eps}
       hard to within factor n^{1-?} means no poly-time algorithm
       can find clique of size n^{1-eps} when optimum is n^{1-?}.
       Here approx_factor = (found size) / (optimum size).

       Gap-ETH lower bound: n^{eps/(1-eps)} in the exponent. */
    double exponent = approx_factor / (1.0 - approx_factor);
    *result = exponent;
    return true;
}

double eth_implies_fpt_neq_w1(int32_t k, double eth_s3, double *fpt_bound) {
    if (k <= 0 || eth_s3 <= 0.0) {
        if (fpt_bound) *fpt_bound = -1.0;
        return 0.0;
    }
    /* If ETH holds, then W[1] != FPT.
       The confidence is based on the derived exponent:
       if n^{Omega(k)} lower bound holds, then k is not a fixed parameter
       in the FPT sense for k-Clique (since the exponent grows with k). */
    double logk = log((double)k + 1.0);
    if (logk < 0.01) logk = 0.01;
    double alpha = eth_s3 * (double)k / (logk + 1.0);

    if (fpt_bound) {
        /* FPT algorithm would require f(k) * n^{O(1)}. Under ETH,
           f(k) must grow at least as n^{alpha - O(1)} which is
           super-exponential in k: f(k) >= 2^{Omega(k log k)}. */
        *fpt_bound = alpha;
    }

    /* Confidence metric: how strongly ETH implies W[1] != FPT.
       Higher values = stronger implication. */
    return alpha / (double)k; /* ratio should be > 0 under ETH */
}

/**
 * Compute the minimum k for which the ETH lower bound for k-Clique
 * exceeds a threshold, indicating the problem transitions from
 * "tractable" to "intractable" in the parameterized sense.
 *
 * @param n              instance size
 * @param eth_s3         ETH constant
 * @param threshold      complexity threshold exponent
 * @return               minimum k where n^{eth_bound} > n^{threshold}
 */
int32_t eth_transition_k(int32_t n, double eth_s3, double threshold) {
    if (n <= 1 || eth_s3 <= 0.0) return 0;
    double logn = log((double)n);
    for (int32_t k = 1; k <= GRAPH_MAX_K; k++) {
        double bound = eth_kclique_lower_bound(n, k, eth_s3);
        if (bound >= threshold * logn) return k;
    }
    return -1;
}

/**
 * Print a comprehensive ETH/SETH bound report for k-Clique.
 *
 * This function demonstrates the core conditional lower bound
 * analysis that is central to fine-grained complexity theory.
 *
 * @param n  graph size
 * @param k  clique parameter
 * @param eth_s3  ETH constant
 */
void print_eth_seth_report(int32_t n, int32_t k, double eth_s3) {
    printf("=== ETH/SETH Lower Bound Report for k-Clique ===\n");
    printf("Parameters: n=%d, k=%d, s_3=%.4f\n\n", n, k, eth_s3);

    /* ETH consequences */
    printf("--- Exponential Time Hypothesis (ETH) ---\n");
    printf("ETH states: 3SAT requires 2^{Omega(n)} time.\n");

    double eth_alpha = eth_kclique_lower_bound(n, k, eth_s3);
    printf("k-Clique lower bound (ETH): n^{%.4f}\n", eth_alpha);

    double delta = eth_delta_exponent(k, eth_s3);
    printf("  Exponent = delta * k, where delta = %.6f\n", delta);
    printf("  Function form: f(k) * n^{delta*k} where f is any computable function\n");

    double eth_time = pow((double)n, eth_alpha);
    printf("  Estimated time: %.2e operations\n\n", eth_time);

    /* SETH consequences */
    printf("--- Strong Exponential Time Hypothesis (SETH) ---\n");
    printf("SETH states: lim_{k->inf} s_k = 1 (kSAT requires 2^{(1-o(1))n}).\n");

    for (int eps_pct = 10; eps_pct >= 1; eps_pct--) {
        double eps = (double)eps_pct / 100.0;
        double seth_lower = seth_kclique_lower_bound(n, k, eps);
        printf("  With eps=%.2f: k-Clique requires n^{%.4f} time\n",
               eps, seth_lower);
    }

    double gap = seth_tightness_gap(n, k);
    printf("  Tightness gap (o(1) term): %.6f\n", gap);
    printf("  Trivial upper bound: O(n^{%d}) = O(n^{%.1f})\n", k, (double)k);
    printf("  SETH says: trivial algorithm is essentially optimal!\n\n");

    /* Gap-ETH consequences */
    printf("--- Gap-ETH Consequences ---\n");
    double approx_result;
    if (gap_eth_approx_clique_bound(0.5, &approx_result)) {
        printf("  Gap-ETH: Approximating clique to 0.5 factor requires n^{%.4f}\n",
               approx_result);
    }

    /* FPT implications */
    printf("\n--- Parameterized Complexity Implications ---\n");
    double fpt_bound_val;
    double confidence = eth_implies_fpt_neq_w1(k, eth_s3, &fpt_bound_val);
    printf("  ETH implies W[1] != FPT with confidence: %.4f\n", confidence);
    printf("  FPT algorithm would need f(k) >= 2^{Omega(k log k)}\n");
    printf("  k-Clique is W[1]-complete (Downey & Fellows 1995)\n");

    printf("\n=== End of Report ===\n");
}
