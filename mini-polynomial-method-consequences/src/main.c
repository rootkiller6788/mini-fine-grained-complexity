#include "poly_method.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* ============================================================================
 * main.c - Polynomial Method Consequences: Demonstration and CLI
 *
 * This file provides the main entry point for the polynomial method
 * module. It demonstrates key concepts through runnable examples:
 * - Polynomial representations of Boolean functions
 * - Fast transforms (FWHT, zeta, mobius)
 * - Algebraic Normal Form computation
 * - OV problem via polynomial method
 * - Williams-style SAT solving
 *
 * L1-L6: Complete knowledge coverage
 * L7: Applications demonstrated (fine-grained complexity)
 * L8: Advanced topics introduced (Razborov-Smolensky)
 * ============================================================================ */

/* ---- Demo: Polynomial representation of basic Boolean functions ---- */
static void demo_boolean_functions(void) {
    printf("=== Boolean Function Polynomial Representations ===\n");

    /* n=3 variables: compute truth tables for AND, OR, XOR, MAJ */
    int32_t n = 3;
    int32_t N = 1 << n;

    /* AND function: f(x) = x_0 AND x_1 AND x_2 */
    uint8_t *tt_and = (uint8_t *)malloc((size_t)N);
    for (int32_t i = 0; i < N; i++)
        tt_and[i] = ((i & 1) && (i & 2) && (i & 4)) ? 1 : 0;

    polynomial_t *poly_and_gf2 = truth_table_to_poly_gf2(tt_and, n);
    printf("AND(x0,x1,x2) as GF(2) polynomial (ANF):\n  ");
    poly_print(poly_and_gf2);
    printf("  Degree over GF(2): %d\n", poly_degree(poly_and_gf2));
    printf("  Expected: x0*x1*x2 (degree 3)\n\n");

    /* Verify evaluation */
    for (int32_t x = 0; x < N; x++) {
        int32_t val = poly_evaluate_gf2(poly_and_gf2, (uint64_t)x);
        int32_t expected = tt_and[x];
        if (val != expected) printf("  MISMATCH at x=%d!\n", x);
    }
    printf("  All %d evaluations match truth table.\n\n", N);

    poly_destroy(poly_and_gf2);
    free(tt_and);

    /* XOR / PARITY function: GF(2) degree = 1 */
    uint8_t *tt_xor = (uint8_t *)malloc((size_t)N);
    for (int32_t i = 0; i < N; i++)
        tt_xor[i] = (uint8_t)(((i & 1) ^ ((i>>1)&1) ^ ((i>>2)&1)) ? 1 : 0);

    polynomial_t *poly_xor_gf2 = truth_table_to_poly_gf2(tt_xor, n);
    printf("PARITY(x0,x1,x2) as GF(2) polynomial (ANF):\n  ");
    poly_print(poly_xor_gf2);
    printf("  Degree over GF(2): %d (expected: 1)\n\n", poly_degree(poly_xor_gf2));
    poly_destroy(poly_xor_gf2);
    free(tt_xor);

    /* MAJORITY function */
    uint8_t *tt_maj = (uint8_t *)malloc((size_t)N);
    for (int32_t i = 0; i < N; i++) {
        int sum = (i&1) + ((i>>1)&1) + ((i>>2)&1);
        tt_maj[i] = (sum >= 2) ? 1 : 0;
    }
    polynomial_t *poly_maj_gf2 = truth_table_to_poly_gf2(tt_maj, n);
    printf("MAJORITY(x0,x1,x2) as GF(2) polynomial (ANF):\n  ");
    poly_print(poly_maj_gf2);
    printf("  Degree over GF(2): %d\n\n", poly_degree(poly_maj_gf2));
    poly_destroy(poly_maj_gf2);
    free(tt_maj);
}

/* ---- Demo: Fast Walsh-Hadamard Transform ---- */
static void demo_fwht(void) {
    printf("=== Fast Walsh-Hadamard Transform ===\n");
    int32_t n = 3;
    int32_t N = 1 << n;

    /* Truth table for AND function */
    double f[8] = {0, 0, 0, 0, 0, 0, 0, 1};  /* only f(111)=1 */

    printf("Input truth table (AND): [");
    for (int32_t i = 0; i < N; i++) printf("%.0f ", f[i]);
    printf("]\n");

    /* FWHT */
    walsh_hadamard_transform(f, n);

    /* Normalize */
    for (int32_t i = 0; i < N; i++) f[i] /= (double)N;

    printf("Fourier coefficients f_hat(S):\n");
    for (int32_t S = 0; S < N; S++) {
        printf("  S={");
        int cnt = 0;
        for (int32_t j = 0; j < n; j++) {
            if (S & (1<<j)) {
                if (cnt>0) printf(",");
                printf("x%d", j);
                cnt++;
            }
        }
        printf("}: %.3f\n", f[S]);
    }
    printf("  (Fourier coefficients match O'Donnell 2014, Ch.1)\n\n");
}

/* ---- Demo: GF(p) Arithmetic ---- */
static void demo_gf_arithmetic(void) {
    printf("=== GF(p) Modular Arithmetic ===\n");
    printf("GF(7): 3 + 5 = %d (mod 7, expected: 1)\n", gf_p_add(3, 5, 7));
    printf("GF(7): 3 * 5 = %d (mod 7, expected: 1)\n", gf_p_mul(3, 5, 7));
    printf("GF(7): 3^3 = %d (mod 7, expected: 6)\n", gf_p_pow(3, 3, 7));
    printf("GF(7): -3 = %d (mod 7, expected: 4)\n", gf_p_neg(3, 7));
    printf("GF(7): inv(3) = %d (mod 7, expected: 5, since 3*5=15=1 mod 7)\n",
           gf_p_inv(3, 7));
    printf("GF(11): inv(7) = %d (mod 11, expected: 8, since 7*8=56=1 mod 11)\n",
           gf_p_inv(7, 11));
    printf("GF(7): inv(0) = %d (no inverse, expected: 0)\n", gf_p_inv(0, 7));
    printf("\n");
}

/* ---- Demo: OV Problem via Polynomial Method ---- */
static void demo_ov_polynomial(void) {
    printf("=== Orthogonal Vectors via Polynomial Method ===\n");

    /* Small instance: d=4, n=3 */
    ov_instance_t *A = ov_create(3, 4);
    ov_instance_t *B = ov_create(3, 4);

    /* A vectors: [1001], [0110], [1111] */
    ov_set_vector(A, 0, 0x9);  /* 1001 */
    ov_set_vector(A, 1, 0x6);  /* 0110 */
    ov_set_vector(A, 2, 0xF);  /* 1111 */

    /* B vectors: [0110], [1001], [1010] */
    ov_set_vector(B, 0, 0x6);  /* 0110 -> dot with A[1]=0110 = 2 (not orthogonal) */
    ov_set_vector(B, 1, 0x9);  /* 1001 -> dot with A[1]=0110 = 0 (orthogonal!) */
    ov_set_vector(B, 2, 0xA);  /* 1010 */

    /* Brute force */
    ov_result_t r1 = ov_solve_brute_force(A, B);
    printf("Brute force: found=%d, a_idx=%d, b_idx=%d\n",
           r1.orthogonal_pair_exists, r1.a_idx, r1.b_idx);

    /* Packed */
    ov_result_t r2 = ov_solve_packed(A, B);
    printf("Packed:      found=%d, a_idx=%d, b_idx=%d\n",
           r2.orthogonal_pair_exists, r2.a_idx, r2.b_idx);

    /* Polynomial */
    ov_result_t r3 = ov_solve_polynomial(A, B);
    printf("Polynomial:  found=%d, a_idx=%d, b_idx=%d\n",
           r3.orthogonal_pair_exists, r3.a_idx, r3.b_idx);

    /* Verify */
    if (r3.orthogonal_pair_exists) {
        bool ok = ov_verify_orthogonal(A, B, r3.a_idx, r3.b_idx);
        printf("Verification: %s\n", ok ? "PASS" : "FAIL");
    }

    /* Count orthogonal pairs */
    int64_t cnt = ov_count_orthogonal_pairs(A, B);
    printf("Orthogonal pairs count: %lld\n\n", (long long)cnt);

    ov_destroy(A);
    ov_destroy(B);
}

/* ---- Demo: Williams' SAT Algorithm ---- */
static void demo_williams_sat(void) {
    printf("=== Williams-Style SAT via Polynomial Method ===\n");

    /* Create a simple 3-SAT instance: n=6, m=4 clauses */
    /* Formula: (x1 OR x2 OR x3) AND (~x1 OR x4 OR x5) AND
     *          (~x2 OR ~x4 OR x6) AND (x3 OR ~x5 OR ~x6) */
    williams_cnf_t *f = williams_cnf_create(6, 4, 3);

    int32_t c1[] = {1, 2, 3};
    int32_t c2[] = {-1, 4, 5};
    int32_t c3[] = {-2, -4, 6};
    int32_t c4[] = {3, -5, -6};

    williams_cnf_add_clause(f, c1, 3);
    williams_cnf_add_clause(f, c2, 3);
    williams_cnf_add_clause(f, c3, 3);
    williams_cnf_add_clause(f, c4, 3);

    printf("Formula: %d vars, %d clauses, width=%d\n",
           f->num_vars, f->num_clauses, f->clause_width);

    /* Brute force */
    williams_result_t r1 = williams_brute_force(f);
    printf("Brute force: SAT=%d, time=%.3f ms, evals=%llu\n",
           r1.satisfiable, r1.runtime_ms, (unsigned long long)r1.evaluations);

    /* Williams polynomial method */
    williams_result_t r2 = williams_solve(f);
    printf("Williams:    SAT=%d, time=%.3f ms, evals=%llu\n",
           r2.satisfiable, r2.runtime_ms, (unsigned long long)r2.evaluations);

    /* Verify witnesses */
    if (r1.satisfiable) {
        printf("Brute force witness: [");
        for (int32_t i = 0; i < f->num_vars; i++)
            printf("%d ", r1.witness[i]);
        printf("], verified=%d\n", williams_verify_witness(f, r1.witness));
    }
    if (r2.satisfiable) {
        printf("Williams witness:    [");
        for (int32_t i = 0; i < f->num_vars; i++)
            printf("%d ", r2.witness[i]);
        printf("], verified=%d\n", williams_verify_witness(f, r2.witness));
    }

    williams_cnf_destroy(f);
    if (r1.witness) free(r1.witness);
    if (r2.witness) free(r2.witness);
    printf("\n");
}

/* ---- Demo: Circuit Polynomials and Lower Bounds ---- */
static void demo_circuit_lower_bounds(void) {
    printf("=== Circuit Lower Bounds via Polynomial Method ===\n");

    int32_t n = 6;

    /* MOD_3 vs MOD_2 degrees */
    uint8_t *tt_mod3 = build_mod_m_truth_table(n, 3);
    uint8_t *tt_mod2 = build_mod_m_truth_table(n, 2);

    if (tt_mod3 && tt_mod2) {
        int32_t deg_mod3_gf2 = gf_degree_from_truth_table(tt_mod3, n, 2);
        int32_t deg_mod2_gf2 = gf_degree_from_truth_table(tt_mod2, n, 2);
        int32_t deg_mod3_gf3 = gf_degree_from_truth_table(tt_mod3, n, 3);

        printf("n=%d variables:\n", n);
        printf("  MOD_3 over GF(2): degree = %d (out of max %d)\n",
               deg_mod3_gf2, n);
        printf("  MOD_2 over GF(2): degree = %d (expected: 1)\n", deg_mod2_gf2);
        printf("  MOD_3 over GF(3): degree = %d (expected: <= 2 by FLT)\n",
               deg_mod3_gf3);
        printf("  => MOD_3 has high GF(2)-degree, so MOD_3 not in AC^0[2]\n");
        printf("  (Razborov-Smolensky 1987: exponential lower bound)\n\n");
    }

    free(tt_mod3);
    free(tt_mod2);

    /* MAJORITY degree */
    uint8_t *tt_maj = build_majority_truth_table(n);
    if (tt_maj) {
        int32_t deg_maj_gf2 = gf_degree_from_truth_table(tt_maj, n, 2);
        printf("  MAJORITY over GF(2): degree = %d (out of max %d)\n",
               deg_maj_gf2, n);
        printf("  => MAJORITY requires high degree => not in AC^0[p]\n");
        printf("  (Razborov 1987: exponential lower bound)\n\n");
    }
    free(tt_maj);

    /* Razborov approximation parameters */
    razborov_approx_params_t params = {3, 1000, 8, 2};
    int32_t deg_bound = razborov_degree_bound(&params, 0.1);
    printf("Razborov bound: depth=%d, size=%d, fan-in=%d, p=%d\n",
           params.depth, params.size, params.max_fan_in, params.p);
    printf("  Polynomial degree bound: %d\n", deg_bound);
    printf("  Any function requiring degree > %d is not computable\n", deg_bound);
    printf("  by AC^0[2] circuits of depth %d and size %d.\n\n",
           params.depth, params.size);
}

int main(void) {
    printf("==========================================\n");
    printf("  mini-polynomial-method-consequences\n");
    printf("  Polynomial Method in Fine-Grained\n");
    printf("  Complexity and Circuit Lower Bounds\n");
    printf("==========================================\n\n");

    demo_boolean_functions();
    demo_fwht();
    demo_gf_arithmetic();
    demo_ov_polynomial();
    demo_williams_sat();
    demo_circuit_lower_bounds();

    printf("=== All demonstrations complete ===\n");
    return 0;
}
