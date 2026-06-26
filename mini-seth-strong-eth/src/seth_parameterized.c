/* seth_parameterized.c -- ETH-based Parameterized Lower Bounds (L7, L8)
 *
 * Connects SETH/ETH to parameterized complexity:
 * - ETH implies FPT != W[1] (k-Clique not FPT)
 * - SETH gives tight lower bounds for FPT algorithms
 * - Connections to treewidth-based algorithms
 * - Lower bounds for dynamic programming on tree decompositions
 */
#include "seth.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

/* ===================================================================
 * ETH ? Parameterized Hardness
 *
 * If ETH holds (3-SAT needs 2^{Omega(n)} time), then:
 * - k-Clique is not in FPT (needs n^{Omega(k)} time)
 * - k-Dominating Set is not in FPT
 * - Many W[1]-hard problems require n^{Omega(k)}
 *
 * This is because the standard parameterized reduction from
 * k-SAT to k-Clique blows up k ? poly(k).
 * =================================================================== */

/* Under ETH, k-Clique parameterized by k requires n^{Omega(k)} time.
 * More precisely: no algorithm with running time f(k) * n^{o(k)}
 * unless ETH is false.
 *
 * Best known: O(n^{omega*k/3}) where omega = matrix mult exponent. */
double eth_kclique_lower_bound(int32_t n, int32_t k) {
    /* ETH ? k-Clique: the reduction from 3-SAT to k-Clique
     * produces an instance with n' = O(2^{n/k}) vertices and
     * parameter k' = O(k).
     * If k-Clique can be solved in f(k) * n'^{g(k)}, then
     * 3-SAT can be solved in 2^{O(g(k)*n/k)} time.
     * For ETH to hold, we need g(k) = Omega(k). */
    if (k <= 0 || n <= 0) return 0.0;
    return pow((double)n, (double)k * 0.1);
}

/* Under ETH, the n^{O(k)} algorithm for k-Dominating Set is optimal.
 * Any n^{o(k)} algorithm would refute ETH. */
double eth_kdomset_lower_bound(int32_t n, int32_t k) {
    if (k <= 0 || n <= 0) return 0.0;
    return pow((double)n, (double)k * 0.1);
}

/* ETH-based lower bound for Hitting Set (k, n).
 * Hitting Set parameterized by solution size k is W[2]-complete.
 * Under ETH: requires n^{Omega(k)} time. */
double eth_hittingset_lower_bound(int32_t n, int32_t k) {
    if (k <= 0 || n <= 0) return 0.0;
    /* W[2]-hardness + ETH => n^{Omega(k)} lower bound */
    return pow((double)n, (double)k * log((double)k) / log(2.0));
}

/* ===================================================================
 * SETH ? Tight FPT Lower Bounds
 *
 * SETH gives more precise lower bounds than ETH for FPT algorithms.
 * For example:
 * - k-Dominating Set needs n^{k-o(1)} under SETH
 * - This is tight: n^{k+o(1)} algorithms exist
 *
 * The n^{k-o(1)} lower bound is EXACT (matching known algorithms),
 * not just asymptotic like ETH gives.
 * =================================================================== */

/* k-Dominating Set: under SETH, needs n^{k-o(1)} time.
 *
 * Reduction: k-SAT ? k-Dominating Set
 * Group n variables into k groups of size n/k.
 * Each group becomes a "choice" gadget.
 * The reduction to Dominating Set with solution size k
 * produces n' = O(k * 2^{n/k}) vertices.
 *
 * If DS can be solved in O(n'^{k-eps}), then k-SAT can be
 * solved in 2^{(1-eps/k)n} time, refuting SETH. */
double seth_kdomset_tight_bound(int32_t n, int32_t k) {
    if (k <= 0 || n <= 0) return 0.0;
    /* SETH gives tight bound: T(n,k) = n^{k-o(1)}
     * The "o(1)" is in the exponent of n as n grows.
     * For a fixed n, the bound is n^{k}. */
    return pow((double)n, (double)k);
}

/* Set Cover with sets of bounded size d:
 * Under SETH, needs (2-eps)^n time for unbounded d.
 * For d=2: O(1.26^n) is possible (polynomial space).
 *
 * These bounds connect SETH to exact exponential algorithms
 * for NP-hard problems beyond SAT. */
double seth_setcover_bound(int32_t n, int32_t d) {
    if (d < 2 || n <= 0) return 0.0;
    /* Set Cover with sets of size ? d:
     * Best known: O((2 - c_d)^n) where c_d ? 0 as d grows.
     * Under SETH: c_d must ? 0.
     * For specific d:
     * d=2 (Vertex Cover): O(1.21^n) best known
     * d=3: O(1.68^n)
     * d=large: approaches 2^n */
    double c_d = 1.0 / pow(2.0, (double)d);
    return pow(2.0 - c_d, (double)n);
}

/* ===================================================================
 * Treewidth and ETH
 *
 * Courcelle's Theorem: MSO properties can be decided in
 * f(tw) * n time on graphs of treewidth ? tw.
 * But f(tw) can be astronomical.
 *
 * Under ETH, f(tw) must be at least 2^{Omega(tw)} for many
 * problems (e.g., 3-Coloring, Vertex Cover on graphs of
 * bounded treewidth).
 * =================================================================== */

/* ETH-optimal treewidth algorithm runtime.
 * For many problems (e.g., Vertex Cover, 3-Coloring, Dominating Set),
 * the optimal dependence on treewidth is:
 *   2^{Theta(tw)} * n^{O(1)}
 *
 * SETH refinement: the base of the exponent matters.
 * e.g., Vertex Cover: O(2^{tw} * n) is optimal under SETH.
 * Dominating Set: O(3^{tw} * n) is optimal under SETH. */
double eth_treewidth_bound(int32_t tw, double base) {
    if (tw < 0) return 0.0;
    return pow(base, (double)tw);
}

/* Specific problem lower bounds under SETH+treewidth:
 *
 * Vertex Cover: base = 2 (matching O(2^{tw} * n) algorithm)
 * Dominating Set: base = 3 (matching O(3^{tw} * n) algorithm)
 * Feedback Vertex Set: base = 3
 * 3-Coloring: base = 3
 * Independent Set: base = 2
 *
 * Each base is tight under SETH. Achieving a smaller base
 * for any of these would refute SETH. */
typedef struct {
    const char *problem_name;
    double      seth_base;    /* optimal base under SETH */
    double      best_known_base;
} treewidth_bound_t;

static const treewidth_bound_t tw_bounds[] = {
    {"Vertex Cover",       2.0, 2.0},
    {"Dominating Set",     3.0, 3.0},
    {"Independent Set",    2.0, 1.9999}, /* 2.0 is optimal under SETH */
    {"Feedback Vertex Set", 3.0, 3.0},
    {"3-Coloring",         3.0, 3.0},
    {"Hamiltonian Cycle",  4.0, 4.0},
    {"Steiner Tree",       3.0, 3.0},
    {"Odd Cycle Transversal", 3.0, 3.0},
    {NULL, 0.0, 0.0}
};

/* Check if a proposed treewidth-based algorithm refutes SETH */
bool check_treewidth_seth_violation(const char *problem,
                                     double proposed_base) {
    for (int32_t i = 0; tw_bounds[i].problem_name != NULL; i++) {
        if (strcmp(tw_bounds[i].problem_name, problem) == 0) {
            /* SETH violation if we beat the SETH-predicted optimal base */
            return proposed_base < tw_bounds[i].seth_base - 0.01;
        }
    }
    return false;
}

/* Print known treewidth-SETH bounds */
void print_treewidth_seth_bounds(void) {
    printf("Treewidth Bounds under SETH:\n");
    printf("%-25s %8s %12s\n", "Problem", "SETH Base", "Best Known");
    printf("----------------------------------------------\n");
    for (int32_t i = 0; tw_bounds[i].problem_name != NULL; i++) {
        printf("%-25s %8.4f %12.4f\n",
               tw_bounds[i].problem_name,
               tw_bounds[i].seth_base,
               tw_bounds[i].best_known_base);
    }
}

/* ===================================================================
 * ETH ? NP-intermediate Problems
 *
 * Under ETH (P != NP implied, but stronger):
 * There exist NP-intermediate problems that are:
 * - Not in P
 * - Not NP-complete (under polynomial reductions)
 *
 * Ladner's theorem constructs such problems unconditionally,
 * but ETH gives more refined information about their complexity.
 *
 * For example: under ETH, certain problems in BQP (quantum
 * polynomial time) might not be in P.
 * =================================================================== */

/* Graph Isomorphism: in NP, not known NP-complete.
 * Under ETH: if GI is NP-complete, then the polynomial
 * hierarchy collapses. Babai (2016): GI in quasipolynomial time.
 *
 * ETH constrains where GI can sit in the complexity landscape. */
const char *eth_graph_isomorphism_status(void) {
    /* Under ETH: GI is NOT NP-complete (unless PH collapses).
     * Babai's algorithm: 2^{O((log n)^3)}.
     * Quasipolynomial but not quite polynomial.
     * ETH would be strained if GI were NP-complete. */
    return "Not NP-complete under ETH (unless PH collapses)";
}

/* Factoring: in NP ? coNP, not known to be in P.
 * Under ETH: unlikely to be NP-complete.
 * Shor's algorithm: O((log n)^3) on quantum computer. */
const char *eth_factoring_status(void) {
    return "Not NP-complete under ETH; BQP but not known in P";
}

/* ===================================================================
 * Computational Complexity Dashboard
 *
 * Collect all SETH/ETH implications into a single report.
 * =================================================================== */

typedef struct {
    const char *problem;
    const char *complexity_class;
    double      lower_bound_exponent;  /* under SETH/ETH */
    const char *algorithm;             /* best known */
    double      upper_bound_exponent;
} problem_complexity_t;

void print_complexity_dashboard(void) {
    printf("\n");
    printf("=== Fine-Grained Complexity Dashboard ===\n");
    printf("%-30s %6s %10s %10s\n", "Problem", "Class", "LB(exp)", "UB(exp)");
    printf("------------------------------------------------------\n");

    /* SAT */
    printf("%-30s %6s %10.4f %10.4f  (PPSZ)\n",
           "3-SAT", "NP-c", 0.0, 0.3863);
    printf("%-30s %6s %10.4f %10.4f  (PPSZ)\n",
           "k-SAT (SETH limit)", "NP-c", 1.0, 1.0);

    /* Graph problems */
    printf("%-30s %6s %10s %10s\n",
           "Vertex Cover", "NP-c", "1.0(O(2^n))", "1.0");
    printf("%-30s %6s %10s %10s\n",
           "k-Clique", "W[1]", "n^{k}", "n^{wk/3}");

    /* String problems under SETH */
    printf("%-30s %6s %10s %10s\n",
           "Edit Distance", "P", "n^{2-o(1)}", "n^2");
    printf("%-30s %6s %10s %10s\n",
           "LCS", "P", "n^{2-o(1)}", "n^2");
    printf("%-30s %6s %10s %10s\n",
           "Frechet Distance", "P", "n^{2-o(1)}", "n^2 log n");

    /* Dynamic problems */
    printf("%-30s %6s %10s %10s\n",
           "APSP", "P", "n^{3-o(1)}", "n^3");
    printf("%-30s %6s %10s %10s\n",
           "3SUM", "P", "n^{2-o(1)}", "n^2");

    printf("------------------------------------------------------\n");
    printf("LB = lower bound under SETH/ETH conjecture\n");
    printf("UB = best known algorithm\n");
}
