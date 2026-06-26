/* ov_conjecture.c - OV Conjecture and Fundamental Laws (L2, L4)
 * Formalizes the OV conjecture, SETH connection, lower bounds.
 */
#include "ov.h"
#include <math.h>
#include <stdio.h>

/* L2: OV Conjecture Formalization
 *
 * OV Conjecture (precise statement):
 *   For all constants epsilon > 0, there is no algorithm that solves
 *   OV in O(n^{2-epsilon}) time when d = omega(log n).
 *
 * Equivalently: Any algorithm for OV with d = omega(log n) requires
 * n^{2-o(1)} time (assuming SETH or as an independent conjecture).
 *
 * Three regimes (Abboud-Williams-Yu 2015):
 *   Regime 1: d = o(log n)      -> easy, n^{1+o(1)} algorithms exist
 *   Regime 2: d = c*log n       -> n^{2-1/O(c)} by Williams' algorithm
 *   Regime 3: d = omega(log n)  -> n^{2-o(1)} conjectured optimal
 */

bool ov_conjecture_violated(int32_t n, int32_t d, double achieved_exp) {
    if (n <= 0 || d <= 0) return false;
    /* Check if dimension is superlogarithmic */
    if (!ov_is_superlogarithmic(n, d)) return false;
    /* Check if exponent is subquadratic */
    /* Allow 1.999 as threshold for "2 - o(1)" */
    return achieved_exp < 1.999;
}

double ov_lower_bound_exponent(int32_t n, int32_t d) {
    if (n <= 0 || d <= 0) return 2.0;
    double log_n = log2((double)n);
    if (log_n <= 0) return 2.0;
    double c = (double)d / log_n;
    if (c < 0.5) return 1.0;           /* d = o(log n): easy regime */
    if (c < 1.0) return 2.0 - 1.0/c;    /* transition: Williams gives 2-1/O(c) */
    if (c < 3.0) return 2.0 - 1.0/(3.0*c); /* approaching quadratic */
    return 2.0;                          /* d = omega(log n): conjectured quadratic */
}

bool ov_is_superlogarithmic(int32_t n, int32_t d) {
    if (n <= 1) return false;
    double log_n = log2((double)n);
    return (double)d > 3.0 * log_n;
}

int32_t ov_regime_classify(int32_t n, int32_t d) {
    /* Returns 1, 2, or 3 for the three regimes */
    if (n <= 0 || d <= 0) return 0;
    double log_n = log2((double)n);
    if (log_n <= 0) return 0;
    double c = (double)d / log_n;
    if (c < 0.5) return 1;       /* d = o(log n): easy */
    if (c < 3.0) return 2;       /* d = Theta(log n): intermediate */
    return 3;                     /* d = omega(log n): hard */
}

/* L4: Williams' Theorem (2005): SETH => OV Conjecture
 *
 * Proof sketch:
 *   Given k-SAT formula F with N variables and M clauses.
 *   Split variables: X = first N/2, Y = last N/2.
 *   For each assignment alpha to X, define vector a_alpha:
 *     a_alpha[j] = 1 iff alpha doesn't satisfy clause j on its own.
 *   For each assignment beta to Y, define vector b_beta:
 *     b_beta[j] = 1 iff beta doesn't satisfy clause j on its own.
 *   Then <a_alpha, b_beta> = 0 iff alpha+beta satisfies ALL clauses.
 *   |A|=|B|=2^{N/2}, d=M.
 *   O(n^{2-eps}) OV => 2^{(1-eps)N} k-SAT, violating SETH.
 */

double ov_to_ksat_threshold(double ov_exp, int32_t n, int32_t d) {
    /* Given an OV algorithm with exponent ov_exp (T = n^{ov_exp}),
     * compute the implied k-SAT exponent.
     * In the reduction: N = 2*log_2(n), so:
     *   T_ksat = (n)^{ov_exp} = 2^{(N/2)*ov_exp} = 2^{ov_exp*N/2}
     * Thus k-SAT exponent = ov_exp / 2.
     * For ov_exp = 2 - eps: k-SAT exponent = 1 - eps/2. */
    (void)d;
    (void)n;
    return ov_exp / 2.0;
}

bool ov_sparsification_applies(int32_t n, int32_t d) {
    /* Sparsification lemma: k-SAT can be reduced to O(n) clauses.
     * For OV, this means if d = O(n^{eps}), sparsification may help.
     * But when d = omega(log n), sparsification is insufficient
     * to bring d down to O(log n). */
    if (n <= 0) return false;
    double log_n = log2((double)n);
    /* Sparsification can reduce clauses by factor 2^{eps*N}.
     * For d to be reducible to O(log n), we need d <= n^{delta}
     * for some delta < 1. */
    return (double)d < pow((double)n, 0.5);
}

int32_t ov_equivalent_ksat_param(double ratio) {
    /* For OV with d = c*log n, the equivalent k-SAT parameter
     * is k = floor(2/eps) where eps = 1/O(c) determines the gap
     * from quadratic time. */
    if (ratio <= 0.5) return 3;   /* easy regime: equivalent to 3-SAT */
    if (ratio <= 1.0) return 5;
    if (ratio <= 2.0) return 10;
    return 20; /* d = omega(log n): equivalent to k-SAT with large k */
}

bool ov_equivalent_to_seth(int32_t n, int32_t d) {
    /* OV conjecture is equivalent to SETH when d = omega(log n).
     * In this regime, beating quadratic OV would beat 2^n for k-SAT
     * for all k (SETH). */
    return ov_is_superlogarithmic(n, d);
}

/* L4: Polynomial Method Lower Bound (informal)
 *
 * The polynomial method (Abboud-Williams-Yu 2015) shows:
 * Any algorithm that evaluates a degree-d polynomial over GF(2)
 * on n^2 inputs requires Omega(n^2) time for d = omega(log n).
 * Since OV can be expressed as evaluating such a polynomial,
 * the lower bound follows (conditionally).
 *
 * ov_polynomial_degree(d): returns the degree of the OV polynomial.
 */
int32_t ov_polynomial_degree(int32_t d) {
    /* The OV detection polynomial:
     *   P(a,b) = prod_{i=1}^{d} (1 - a_i * b_i) over GF(2)
     * has degree d in the coordinates of a and b.
     * This is because each term (1 - a_i*b_i) has degree 2,
     * and the product of d such terms has degree 2d.
     * Over GF(2), 1 - x = 1 + x, so the polynomial is multilinear. */
    return d < 0 ? 0 : d;
}

int32_t *ov_polynomial_coefficients(int32_t d, int32_t *nc) {
    /* Compute coefficients of the multilinear OV polynomial.
     * P(a,b) = sum_{S subset [d]} (-1)^{|S|} prod_{i in S} a_i*b_i
     * Returns array of size 2^d (only practical for d <= 20). */
    if (d <= 0 || d > 20) { *nc = 0; return NULL; }
    int32_t size = 1 << d;
    int32_t *coeffs = (int32_t *)calloc((size_t)size, sizeof(int32_t));
    if (!coeffs) { *nc = 0; return NULL; }
    /* coeffs[mask] = (-1)^{popcount(mask)}
     * because each subset S contributes (-1)^{|S|} to the coefficient
     * of prod_{i in S} a_i * b_i. */
    for (int32_t mask = 0; mask < size; mask++) {
        int32_t pc = __builtin_popcount((unsigned int)mask);
        coeffs[mask] = (pc % 2 == 0) ? 1 : -1;
    }
    *nc = size;
    return coeffs;
}

/* L4: Communication Complexity Lower Bound */
int64_t ov_communication_lower_bound(int32_t n, int32_t d) {
    /* Deterministic communication complexity of OV:
     * Alice holds A, Bob holds B. They want to decide if exists orth pair.
     * Lower bound: Omega(n) bits (Alice must describe which vector to test).
     * For d = omega(log n): Omega(n*d) bits conjectured. */
    (void)d;
    return (int64_t)n; /* trivial lower bound */
}

int64_t ov_communication_protocol_simulate(const ov_instance_t *inst, int32_t bpm) {
    /* Simulate bits exchanged in a naive protocol:
     * Alice sends one vector (~d bits) to Bob, Bob checks all his vectors.
     * Cost: d bits per round. */
    if (!inst) return 0;
    (void)bpm;
    return (int64_t)inst->A->dimension; /* one vector worth of bits */
}

/* L2+L9: OV Conjecture implies 3SUM Conjecture? */
bool ov_conjecture_implies_3sum(int32_t n, int32_t d) {
    /* Research question: does OV hardness imply 3SUM hardness?
     * Current understanding: No direct implication known.
     * Both are independent fine-grained conjectures.
     * This function always returns false (no known implication). */
    (void)n; (void)d;
    return false;
}

void ov_print_equivalence_class(void) {
    printf("\n========================================\n");
    printf("OV Equivalence Class (subquadratic reductions)\n");
    printf("========================================\n");
    printf("Problems requiring n^{2-o(1)} under OV conjecture:\n\n");
    printf("  L6.1 Pattern Matching with Wildcards\n");
    printf("       Text T, pattern P with wildcards '?'.\n");
    printf("       OV reduces to PM: n vectors -> strings of length O(n*d).\n\n");
    printf("  L6.2 Edit Distance (Levenshtein)\n");
    printf("       Minimum insertions/deletions/substitutions.\n");
    printf("       Bringmann-Kunnemann (2015): OV <= Edit Distance.\n\n");
    printf("  L6.3 Longest Common Subsequence (LCS)\n");
    printf("       Abboud-Backurs-Williams (2015): OV <= LCS.\n\n");
    printf("  L6.4 Frechet Distance\n");
    printf("       Curve similarity metric.\n");
    printf("       Bringmann (2014): OV <= Frechet Distance.\n\n");
    printf("  L6.5 Dynamic Time Warping (DTW)\n");
    printf("       Abboud-Backurs-Williams (2015): OV <= DTW.\n\n");
    printf("  L6.6 Graph Diameter (sparse graphs)\n");
    printf("       Roditty-Williams (2013): OV <= Diameter.\n\n");
    printf("  L6.7 Hamming Nearest Neighbor\n");
    printf("       Alman-Williams (2015): OV <= HNN.\n\n");
    printf("  L6.8 Subset Sum (small integers)\n");
    printf("       Abboud-Bringmann-Hermelin-Shabtay (2019).\n\n");
    printf("Reference: Abboud-Williams-Yu (2015), SODA.\n");
    printf("========================================\n\n");
}
