#include "poly_method.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ============================================================================
 * lower_bounds.c - Lower Bounds via Polynomial Representations
 *
 * The polynomial method is a powerful tool for proving circuit lower bounds.
 * The core idea: if a Boolean function can be computed by a "small" circuit
 * (in some class), then it can be approximated by a "low-degree" polynomial.
 * Conversely, if a function requires "high-degree" polynomials, it cannot
 * have small circuits.
 *
 * Key results proved via the polynomial method:
 *
 * 1. Razborov (1987): MAJORITY not in AC^0[p] for any prime p.
 *    - MAJORITY requires exponential-size constant-depth circuits
 *      even with MOD_p gates.
 *
 * 2. Smolensky (1987): MOD_q not in AC^0[p] for distinct primes p, q.
 *    - MOD_3 cannot be computed by polynomial-size constant-depth
 *      circuits with MOD_2 (PARITY) gates.
 *
 * 3. Williams (2014): NEXP not in ACC^0.
 *    - Nondeterministic exponential time is not contained in
 *      constant-depth circuits with MOD_m gates. Proved by combining
 *      the polynomial method with algorithm design.
 *
 * L4: Fundamental laws - Razborov-Smolensky lower bounds
 * L5: Algorithms - degree lower bound computation
 * L8: Advanced topics - circuit lower bounds, natural proofs barrier
 * L9: Research frontiers - current barriers, meta-complexity
 * ============================================================================ */

/* ============================================================================
 * Compute the GF(p) Degree of a Boolean Function
 *
 * The GF(p)-degree of f: {0,1}^n -> {0,1} is the minimum degree of
 * a multi-linear polynomial over GF(p) that agrees with f on all
 * Boolean inputs.
 *
 * Theorem: Every Boolean function has a UNIQUE multi-linear polynomial
 * over GF(p) representing it. The degree of this polynomial is the
 * GF(p)-degree of the function.
 *
 * Method: Convert truth table to polynomial over GF(p), then find
 * the maximum monomial degree with non-zero coefficient.
 *
 * Time: O(n * 2^n).
 * ============================================================================ */

int32_t gf_degree_from_truth_table(const uint8_t *truth_table, int32_t n,
                                     int32_t p) {
    if (!truth_table || n < 0 || n > 20) return -1;
    int32_t N = 1 << n;

    if (p == 2) {
        polynomial_t *poly = truth_table_to_poly_gf2(truth_table, n);
        if (!poly) return -1;
        int32_t deg = poly_degree(poly);
        poly_destroy(poly);
        return deg;
    }

    /* For general p: compute the unique multi-linear polynomial over GF(p).
     * Use the Mobius inversion formula over GF(p):
     * f(x) = sum_{S} c_S * prod_{i in S} x_i
     * c_S = sum_{T subseteq S} (-1)^{|S|-|T|} * f(T)  (mod p)
     */
    int32_t *coeffs = (int32_t *)calloc((size_t)N, sizeof(int32_t));
    if (!coeffs) return -1;

    for (int32_t i = 0; i < N; i++)
        coeffs[i] = (int32_t)(truth_table[i] & 1);

    /* Mobius transform over integers */
    for (int32_t i = 0; i < n; i++) {
        for (int32_t mask = 0; mask < N; mask++) {
            if (mask & (1 << i)) {
                coeffs[mask] = (coeffs[mask] - coeffs[mask ^ (1 << i)]) % p;
                if (coeffs[mask] < 0) coeffs[mask] += p;
            }
        }
    }

    /* Find maximum degree with non-zero coefficient */
    int32_t max_deg = 0;
    for (int32_t S = 0; S < N; S++) {
        if (coeffs[S] != 0) {
            int32_t deg = 0;
            for (int32_t i = 0; i < n; i++)
                if (S & (1 << i)) deg++;
            if (deg > max_deg) max_deg = deg;
        }
    }

    free(coeffs);
    return max_deg;
}

/* ============================================================================
 * The MOD_m Function (0/1 valued)
 *
 * MOD_m(x_1,...,x_n) = 1 if sum_i x_i ? 0 (mod m), else 0.
 *
 * Property: Over GF(p) for prime p NOT dividing m, MOD_m requires
 * polynomial degree Omega(n). This is the foundation of the
 * Razborov-Smolensky separation.
 *
 * Special cases:
 * - MOD_2 = PARITY
 * - MOD_2 over GF(2) has degree 1 (just x_1 + x_2 + ... + x_n)
 * - MOD_3 over GF(2) requires degree Omega(n)
 * ============================================================================ */

uint8_t *build_mod_m_truth_table(int32_t n, int32_t m) {
    int32_t N = 1 << n;
    uint8_t *tt = (uint8_t *)malloc((size_t)N * sizeof(uint8_t));
    if (!tt) return NULL;
    for (int32_t x = 0; x < N; x++) {
        int32_t sum = 0;
        for (int32_t i = 0; i < n; i++)
            if (x & (1 << i)) sum++;
        tt[x] = (uint8_t)((sum % m == 0) ? 1 : 0);
    }
    return tt;
}

/* ============================================================================
 * Verify Razborov-Smolensky empirically: MOD_3 over GF(2)
 *
 * Theorem (Smolensky 1987):
 *   Any polynomial over GF(2) that agrees with MOD_3 on {0,1}^n
 *   must have degree Omega(n).
 *
 * We can verify this for small n by computing the exact degree.
 * For n = 8, MOD_3 over GF(2) requires degree >= n/2, while
 * MOD_2 over GF(2) requires degree exactly 1.
 * ============================================================================ */

/* Compute empirical degree of MOD_m over GF(p) for small n */
int32_t mod_m_degree_over_p(int32_t n, int32_t m, int32_t p) {
    uint8_t *tt = build_mod_m_truth_table(n, m);
    if (!tt) return -1;
    int32_t deg = gf_degree_from_truth_table(tt, n, p);
    free(tt);
    return deg;
}

/* ============================================================================
 * The MAJORITY Function
 *
 * MAJ(x_1,...,x_n) = 1 iff sum_i x_i >= n/2.
 *
 * Razborov (1987): MAJORITY is not in AC^0[p] for any prime p.
 * Specifically, any depth-d AC^0[p] circuit for MAJORITY requires
 * size exp(Omega(n^{1/(2d-2)})).
 *
 * The proof uses the polynomial method:
 * 1. Every AC^0[p] circuit of size s and depth d can be approximated
 *    by a polynomial of degree (log s)^{O(d)} over GF(p).
 * 2. MAJORITY cannot be approximated by low-degree polynomials
 *    over any field: any degree-d polynomial agrees with MAJORITY
 *    on at most 1/2 + O(d/sqrt(n)) fraction of inputs.
 * 3. Therefore, s must be exp(Omega(n^{1/O(d)})).
 * ============================================================================ */

uint8_t *build_majority_truth_table(int32_t n) {
    int32_t N = 1 << n;
    uint8_t *tt = (uint8_t *)malloc((size_t)N * sizeof(uint8_t));
    if (!tt) return NULL;
    for (int32_t x = 0; x < N; x++) {
        int32_t sum = 0;
        for (int32_t i = 0; i < n; i++)
            if (x & (1 << i)) sum++;
        tt[x] = (uint8_t)((sum >= (n + 1) / 2) ? 1 : 0);
    }
    return tt;
}

/* ============================================================================
 * Degree Lower Bound for MAJORITY
 *
 * Theorem: Any real polynomial p of degree d that tries to approximate
 * MAJORITY must satisfy:
 *   max_{x in {0,1}^n} |p(x) - MAJ(x)| >= c * d / sqrt(n)
 * for some constant c > 0.
 *
 * In particular, exact representation requires degree n.
 * ============================================================================ */

int32_t majority_degree_over_p(int32_t n, int32_t p) {
    uint8_t *tt = build_majority_truth_table(n);
    if (!tt) return -1;
    int32_t deg = gf_degree_from_truth_table(tt, n, p);
    free(tt);
    return deg;
}

/* ============================================================================
 * The AND-OR Tree / Decision Tree Complexity
 *
 * A decision tree for f: {0,1}^n -> {0,1} is a binary tree where
 * internal nodes query variables and leaves output 0 or 1.
 *
 * Theorem: deg_{GF(2)}(f) <= D(f) (decision tree depth)
 * because the polynomial representation can be read off the
 * decision tree: each leaf corresponds to a monomial.
 *
 * Therefore, functions with high GF(2)-degree require deep
 * decision trees (and hence large circuit depth in certain models).
 * ============================================================================ */

/* Construct the unique GF(2) polynomial for given truth table
 * and return its degree */
int32_t decision_tree_degree_lower_bound(const uint8_t *truth_table,
                                            int32_t n) {
    return gf_degree_from_truth_table(truth_table, n, 2);
}

/* ============================================================================
 * The Polynomial Method and the Natural Proofs Barrier
 *
 * Razborov and Rudich (1997) proved that most known circuit lower bound
 * techniques fit into the "natural proofs" framework and cannot prove
 * strong lower bounds (like P != NP) under standard cryptographic
 * assumptions.
 *
 * The polynomial method, however, is one of the few techniques that
 * CAN sometimes bypass the natural proofs barrier because:
 * 1. It works with algebraic structures (polynomials over finite fields)
 * 2. It leverages properties of specific functions (like MOD_m)
 * 3. It combines circuit analysis with algorithm design
 *
 * Williams' ACC lower bound (NEXP not in ACC^0) is a prime example
 * where the polynomial method succeeded where natural proofs could not.
 *
 * L8: Advanced topics - natural proofs barrier
 * L9: Research frontiers - current lower bound techniques
 * ============================================================================ */

/* Check if a given truth table corresponds to a "natural" property.
 * A natural property (against a circuit class C) is:
 * - Constructive: decidable in time 2^{O(n)}
 * - Largeness: holds for at least a 1/n^{O(1)} fraction of functions
 * - Usefulness: no function in C has the property
 *
 * The natural proofs barrier says: under the existence of exponentially
 * hard pseudo-random generators, no natural property can be used to
 * prove super-polynomial lower bounds against any circuit class
 * capable of computing PRGs.
 *
 * Our check: given a truth table, is it "large" (many functions share
 * the property) but "useful" (easy to compute functions don't have it)?
 * This is a simplified educational check. */

/* Analyze a truth table for natural property characteristics.
 * We use the low-degree criterion: functions with degree <= d
 * are "easy" (computable by small AC^0 circuits, by Razborov's theorem).
 * Functions with degree > d could be "hard."
 * 
 * A property "degree > d" is:
 * - Large if d is small enough (most functions have high degree)
 * - Useful for AC^0 if d = polylog(n) (AC^0 functions have low degree)
 * 
 * The issue: such a property is also constructive (degree computable
 * in O(n*2^n) time), making it a NATURAL property. By the natural
 * proofs barrier, it cannot prove P != NP (but CAN prove lower bounds
 * for weaker classes like AC^0[2], which is exactly what Razborov-
 * Smolensky does). */
natural_property_check_t analyze_truth_table_property(
    const uint8_t *truth_table, int32_t n) {
    natural_property_check_t result = {false, false, 0.0, 0};

    int32_t deg = decision_tree_degree_lower_bound(truth_table, n);

    /* Degree distribution among all Boolean functions:
     * The number of functions with degree <= d over GF(2) is
     * 2^{sum_{i=0}^{d} C(n,i)}, out of 2^{2^n} total.
     * For large n, if d << n/2, this is a tiny fraction. */
    int32_t threshold = n / 3;  /* illustrative */
    result.degree_threshold = threshold;

    result.is_useful = (deg > threshold);
    /* largeness: fraction of functions with degree > threshold */
    /* sum_{i=0}^{threshold} C(n,i) approx */
    int64_t comb_sum = 1;  /* i=0 */
    int64_t comb_nk = 1;
    for (int32_t k = 1; k <= threshold && k <= n; k++) {
        comb_nk = comb_nk * (n - k + 1) / k;
        comb_sum += comb_nk;
    }
    /* Functions with degree <= threshold: 2^{comb_sum} out of 2^{2^n} */
    result.largeness = 1.0 - pow(2.0, (double)(comb_sum - (1 << n)));
    if (result.largeness < 0.0) result.largeness = 0.0;
    result.is_large = (result.largeness > 0.01);  /* at least 1% */

    return result;
}

/* ============================================================================
 * Razborov's Approximation Method
 *
 * Razborov (1987) introduced the method of approximating Boolean circuits
 * by low-degree polynomials over finite fields. The key steps:
 *
 * 1. Replace each gate with a "probabilistic polynomial" that
 *    approximates the Boolean function of the gate.
 * 2. AND gate with fan-in k: approximated by a random degree-D
 *    polynomial over GF(p), error ? k/p^D.
 * 3. OR gate: by De Morgan, OR(x_1,...,x_k) = NOT(AND(NOT(x_i))).
 * 4. After composing d levels, the total degree is D^d and
 *    the total error is at most (total number of AND gates) * k/p^D.
 * 5. Choose D = O(log s) so that degree = poly(log s) and error < 1/3.
 *
 * This function simulates the approximation process for a given
 * circuit structure and computes the resulting polynomial degree bound.
 *
 * L5: Algorithms - Razborov's approximation method
 * L4: Fundamental laws - degree/size tradeoff theorem
 * ============================================================================ */

/* Compute the polynomial degree needed to approximate a circuit
 * with given parameters, achieving error < epsilon. */
int32_t razborov_degree_bound(const razborov_approx_params_t *params,
                                double epsilon) {
    if (!params) return -1;

    int32_t d = params->depth;
    int32_t s = params->size;
    int32_t k = params->max_fan_in;
    int32_t p = params->p;

    /* We need: s * k / p^D < epsilon, degree <= D^d
     * So D > log_p(s * k / epsilon)
     * degree = D^d */

    double min_D = log((double)(s * k) / epsilon) / log((double)p);
    int32_t D = (int32_t)ceil(min_D);
    if (D < 1) D = 1;

    /* Degree bound: D^d */
    int32_t degree = 1;
    for (int32_t i = 0; i < d; i++) {
        if (degree > 1000000) return 1000000;  /* cap for practicality */
        degree *= D;
    }

    return degree;
}

/* ============================================================================
 * Lower Bound Implication: if a function requires degree > B,
 * it cannot be computed by any circuit with the given depth/size
 * parameters. This is the contrapositive of the approximation theorem.
 * ============================================================================ */

bool razborov_lower_bound_check(int32_t function_degree,
                                  const razborov_approx_params_t *params,
                                  double epsilon) {
    int32_t bound = razborov_degree_bound(params, epsilon);
    return function_degree > bound;
}

/* ============================================================================
 * Smolensky's Lower Bound: MOD_q not in AC^0[p]
 *
 * Theorem (Smolensky 1987): For distinct primes p and q, the function
 * MOD_q (which is 1 iff sum of inputs is divisible by q) cannot be
 * computed by AC^0[p] circuits of polynomial size.
 *
 * Proof sketch:
 * 1. Every AC^0[p] circuit of size poly(n) and depth O(1) can be
 *    approximated (error < 1/3) by a degree polylog(n) polynomial
 *    over GF(p).
 * 2. MOD_q over GF(p) requires degree Omega(n) for exact computation.
 * 3. Moreover, MOD_q is "correlated with itself": no low-degree
 *    polynomial over GF(p) can approximate it within error < 1/3.
 * 4. Therefore, MOD_q requires super-polynomial size AC^0[p] circuits.
 *
 * This function demonstrates the gap by computing the exact GF(p)-degree
 * of MOD_q for small n and showing it grows with n.
 * ============================================================================ */

smolensky_result_t smolensky_analyze(int32_t n, int32_t mod_m, int32_t field_p) {
    smolensky_result_t result = {n, mod_m, field_p, 0, false};

    uint8_t *tt = build_mod_m_truth_table(n, mod_m);
    if (!tt) return result;

    result.exact_degree = gf_degree_from_truth_table(tt, n, field_p);
    free(tt);

    /* Heuristic: degree > n/3 indicates "requires high degree" */
    result.requires_high_degree = (result.exact_degree > n / 3);

    return result;
}

/* ============================================================================
 * The Polynomial Method Heuristic for P vs NP
 *
 * While the polynomial method has been successful for circuit lower
 * bounds (AC^0, AC^0[p], ACC^0), extending it to prove P != NP
 * faces fundamental barriers:
 *
 * 1. The natural proofs barrier (Razborov-Rudich 1997): any "natural"
 *    combinatorial property cannot prove strong circuit lower bounds
 *    under cryptographic assumptions.
 *
 * 2. The relativization barrier (Baker-Gill-Solovay 1975): any proof
 *    technique that relativizes cannot resolve P vs NP.
 *
 * 3. The algebrization barrier (Aaronson-Wigderson 2009): extensions
 *    of the polynomial method that "algebrize" still cannot resolve
 *    P vs NP.
 *
 * Current frontier: Meta-complexity (the complexity of computing
 * complexity measures like MCSP) may offer a path forward.
 *
 * L9: Research frontiers - meta-complexity, barriers in complexity theory
 * ============================================================================ */

/* Compute the circuit complexity of a truth table (minimum circuit size).
 * This is a META-COMPUTATIONAL problem: we're computing a property
 * of a Boolean function given its truth table.
 *
 * The Minimum Circuit Size Problem (MCSP):
 *   Input: Truth table of f and a size bound s
 *   Question: Does f have a circuit of size <= s?
 *
 * MCSP is in NP but not known to be NP-complete, and resolving its
 * status is a major open problem. It is believed to be a candidate
 * NP-intermediate problem (by Ladner's theorem, if P != NP).
 *
 * This function estimates circuit complexity by analyzing polynomial degree. */
circuit_complexity_estimate_t estimate_circuit_complexity(
    const uint8_t *truth_table, int32_t n) {
    circuit_complexity_estimate_t est = {n, 1 << n, 0, 0, 0.0};
    if (!truth_table || n < 0) return est;

    polynomial_t *poly = truth_table_to_poly_gf2(truth_table, n);
    if (!poly) return est;

    est.exact_degree = poly_degree(poly);
    est.sparsity = poly_sparsity(poly);

    /* Heuristic circuit size estimate based on polynomial complexity.
     * A circuit can compute each monomial with degree-1 gates,
     * and combine them. The number of gates is roughly proportional
     * to the sparsity times the degree. */
    est.complexity_estimate = (double)est.sparsity * (double)(est.exact_degree + 1);

    poly_destroy(poly);
    return est;
}

/* ============================================================================
 * Polynomial Method and Derandomization
 *
 * The polynomial method also connects to derandomization:
 * - If we can construct a polynomial that distinguishes random from
 *   pseudorandom distributions, we get hardness vs randomness tradeoffs.
 * - Nisan-Wigderson (1994): Hardness vs Randomness framework uses
 *   low-degree polynomial extensions of Boolean functions.
 * - The permanent's random self-reducibility (Lipton 1991) uses
 *   univariate polynomial interpolation.
 *
 * This function checks whether a given truth table can be extended
 * to a low-degree multi-linear polynomial over a larger field,
 * a technique used in interactive proofs (IP=PSPACE, GKR protocol).
 * ============================================================================ */

/* Check if a Boolean function given by its truth table is the restriction
 * to {0,1}^n of a low-degree multi-linear polynomial over GF(p).
 * This is always true (by interpolation), but the degree may be high.
 * Returns the minimum degree needed for such an extension. */
int32_t boolean_extension_degree(const uint8_t *truth_table, int32_t n,
                                    int32_t p) {
    return gf_degree_from_truth_table(truth_table, n, p);
}

/* ============================================================================
 * Razborov-Smolensky: Construct an Explicit High-Degree Function
 *
 * We construct the MOD_3 function on n variables and compute its
 * GF(2)-degree. For n = 8, this is typically 8 (max possible),
 * demonstrating the Omega(n) lower bound phenomenon.
 *
 * Contrast: MOD_2 (PARITY) over GF(2) has degree exactly 1.
 * ============================================================================ */

degree_comparison_t compare_mod_degrees(int32_t n) {
    degree_comparison_t comp = {n, 0, 2, 0, 0.0};

    /* MOD_3 over GF(2) */
    uint8_t *tt3 = build_mod_m_truth_table(n, 3);
    if (tt3) {
        comp.degree = gf_degree_from_truth_table(tt3, n, 2);
        comp.mod_m = 3;
        comp.ratio = (n > 0) ? (double)comp.degree / (double)n : 0.0;
        free(tt3);
    }

    return comp;
}
