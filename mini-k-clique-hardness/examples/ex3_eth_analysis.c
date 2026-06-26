/**
 * ex3_eth_analysis.c ? ETH/SETH Lower Bound Analysis
 *
 * Demonstrates the conditional lower bound methodology:
 * computes k-Clique time lower bounds under ETH and SETH.
 *
 * Shows how the parameter k affects the lower bound exponent
 * and demonstrates the fine-grained complexity analysis.
 *
 * Knowledge: L2/L4 ? ETH/SETH consequences, conditional lower bounds
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/kclique_types.h"
#include "../include/kclique_eth.h"

int main(void) {
    printf("=== ETH/SETH Lower Bound Analysis for k-Clique ===\n\n");

    /* ETH Analysis */
    printf("--- Exponential Time Hypothesis (ETH) ---\n");
    printf("ETH: 3SAT requires 2^{Omega(n)} time (s_3 > 0).\n");
    printf("Consequence: k-Clique requires n^{Omega(k)} time.\n\n");

    double s_3 = 1.0; /* conservative ETH constant */

    printf("ETH k-Clique Lower Bounds (n=1000, s_3=%.1f):\n", s_3);
    printf("  k  | ETH Bound (exponent) | Implied Runtime\n");
    printf("-----+-----------------------+----------------\n");
    for (int32_t k = 2; k <= 10; k++) {
        double alpha = eth_kclique_lower_bound(1000, k, s_3);
        double runtime = pow(1000.0, alpha);
        printf("  %2d | n^{%.4f}          | ~%.2e ops\n",
               k, alpha, runtime);
    }
    printf("\n");

    /* SETH Analysis */
    printf("--- Strong Exponential Time Hypothesis (SETH) ---\n");
    printf("SETH: lim_{k->inf} s_k = 1 (kSAT requires 2^{(1-o(1))n}).\n");
    printf("Consequence: k-Clique requires n^{k-o(1)} time (TIGHT).\n\n");

    printf("SETH k-Clique Lower Bounds (n=1000):\n");
    printf("  k  | SETH Bound (exponent) | Trivial Upper | Gap\n");
    printf("-----+-----------------------+---------------+-----\n");
    for (int32_t k = 2; k <= 10; k++) {
        double seth_lower = seth_kclique_lower_bound(1000, k, 0.01);
        double trivial = (double)k;
        double gap = seth_tightness_gap(1000, k);
        printf("  %2d | n^{%.4f}          | n^{%.1f}         | %.4f\n",
               k, seth_lower, trivial, gap);
    }
    printf("\n");

    /* ETH vs SETH comparison */
    printf("--- ETH vs SETH Comparison ---\n");
    printf("ETH:  k-Clique requires n^{Omega(k)} time\n");
    printf("      (bound grows with k but with a factor < 1)\n");
    printf("      Example: n^{delta*k} for delta = s_3/poly(k)\n\n");
    printf("SETH: k-Clique requires n^{k-o(1)} time\n");
    printf("      (bound is essentially k, matching trivial O(n^k))\n");
    printf("      SETH gives the tightest possible lower bound!\n\n");

    /* Conditional Lower Bounds */
    printf("--- Conditional Lower Bound Catalog ---\n");
    conditional_lower_bound_t clb_eth = make_conditional_lower_bound(
        "ETH", "k-Clique", 2.5, false, false);
    conditional_lower_bound_t clb_seth = make_conditional_lower_bound(
        "SETH", "k-Clique", 4.0, false, true);

    printf("Bound 1: "); print_clb(&clb_eth);
    printf("Bound 2: "); print_clb(&clb_seth);

    int cmp = compare_clb_hardness(&clb_eth, &clb_seth);
    printf("Comparison: SETH bound is ");
    if (cmp < 0) printf("STRONGER (higher exponent)\n");
    else if (cmp > 0) printf("WEAKER\n");
    else printf("EQUIVALENT\n");

    free_clb(&clb_eth);
    free_clb(&clb_seth);

    /* Parameterized complexity implication */
    printf("\n--- Parameterized Complexity Implication ---\n");
    double fpt_bound;
    double confidence = eth_implies_fpt_neq_w1(5, s_3, &fpt_bound);
    printf("Under ETH (s_3=%.1f, k=5):\n", s_3);
    printf("  Confidence that FPT != W[1]: %.4f\n", confidence);
    printf("  FPT algorithm would need f(k) >= 2^{Omega(%.1f)}\n", fpt_bound);

    printf("\nDone.\n");
    return 0;
}
