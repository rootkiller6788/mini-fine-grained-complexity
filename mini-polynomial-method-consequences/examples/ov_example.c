#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "poly_method.h"

/* ============================================================================
 * ov_example.c - Orthogonal Vectors via Polynomial Method
 *
 * L6: Canonical problem - OV detection and counting
 * L7: Application - fine-grained complexity conditional lower bounds
 * ============================================================================ */

int main(void) {
    printf("=== Orthogonal Vectors: Polynomial Method Example ===\n\n");

    int32_t d = 8, n = 6;

    /* Create random OV instances for demonstration */
    ov_instance_t *A = ov_create(n, d);
    ov_instance_t *B = ov_create(n, d);

    uint64_t state = 12345;
    uint64_t mask = (1ULL << d) - 1;
    for (int32_t i = 0; i < n; i++) {
        state = state * 6364136223846793005ULL + 1442695040888963407ULL;
        A->vectors[i] = state & mask;
        state = state * 6364136223846793005ULL + 1442695040888963407ULL;
        B->vectors[i] = state & mask;
    }

    printf("Instance: n=%d, d=%d\n", n, d);
    printf("Set A vectors:\n");
    for (int32_t i = 0; i < n; i++) {
        printf("  A[%d] = ", i);
        for (int32_t j = 0; j < d; j++)
            printf("%d", (int)((A->vectors[i] >> j) & 1));
        printf("\n");
    }
    printf("Set B vectors:\n");
    for (int32_t i = 0; i < n; i++) {
        printf("  B[%d] = ", i);
        for (int32_t j = 0; j < d; j++)
            printf("%d", (int)((B->vectors[i] >> j) & 1));
        printf("\n");
    }

    printf("\n--- Brute Force (O(n^2 * d)) ---\n");
    ov_result_t r1 = ov_solve_brute_force(A, B);
    printf("Orthogonal pair exists: %s\n", r1.orthogonal_pair_exists ? "YES" : "NO");
    if (r1.orthogonal_pair_exists)
        printf("  Pair: A[%d] . B[%d] = 0\n", r1.a_idx, r1.b_idx);

    printf("\n--- Packed (O(n^2), d <= 64) ---\n");
    ov_result_t r2 = ov_solve_packed(A, B);
    printf("Orthogonal pair exists: %s\n", r2.orthogonal_pair_exists ? "YES" : "NO");
    if (r2.orthogonal_pair_exists)
        printf("  Pair: A[%d] . B[%d] = 0\n", r2.a_idx, r2.b_idx);

    /* Verify */
    if (r2.orthogonal_pair_exists) {
        bool ok = ov_verify_orthogonal(A, B, r2.a_idx, r2.b_idx);
        printf("  Verification: %s\n", ok ? "PASS" : "FAIL");
    }

    printf("\n--- Polynomial Method ---\n");
    ov_result_t r3 = ov_solve_polynomial(A, B);
    printf("Orthogonal pair exists: %s\n", r3.orthogonal_pair_exists ? "YES" : "NO");
    if (r3.orthogonal_pair_exists)
        printf("  Pair: A[%d] . B[%d] = 0\n", r3.a_idx, r3.b_idx);

    printf("\n--- Statistics ---\n");
    int64_t total_pairs = ov_count_orthogonal_pairs(A, B);
    printf("Total orthogonal pairs: %lld out of %d possible\n",
           (long long)total_pairs, n * n);
    printf("Density: %.2f%%\n", 100.0 * (double)total_pairs / (double)(n * n));

    /* Demonstrate polynomial construction for a specific vector */
    printf("\n--- Polynomial for B[0] ---\n");
    polynomial_t *pb = poly_create(d, 64, true);
    /* P_b(x) = prod_{i: b_i=1} (1 + x_i) */
    poly_add_term(pb, 1.0, 0ULL);  /* constant 1 */
    uint64_t b0 = B->vectors[0];
    uint64_t sub = b0;
    while (sub) {
        poly_add_term(pb, 1.0, sub);
        sub = (sub - 1) & b0;
    }
    printf("P_b polynomial: degree=%d, terms=%d\n",
           poly_degree(pb), pb->num_terms);
    printf("(P_b has 2^{|b|} terms, where |b| = popcount of B[0])\n");

    /* Check: P_b(a) = 1 iff a . b = 0 */
    printf("\nEvaluation of P_b on all A vectors:\n");
    for (int32_t i = 0; i < n; i++) {
        int32_t val = poly_evaluate_gf2(pb, A->vectors[i]);
        int32_t dot = 0;
        uint64_t prod_val = A->vectors[i] & B->vectors[0];
        for (int32_t j = 0; j < d; j++)
            if (prod_val & (1ULL << j)) dot++;
        printf("  A[%d]: P_b=%d, dot=%d %s\n", i, val, dot,
               (val == 1 && dot == 0) ? "(correct)" :
               (val == 0 && dot > 0) ? "(correct)" : "MISMATCH!");
    }

    poly_destroy(pb);
    ov_destroy(A);
    ov_destroy(B);

    printf("\n=== Example complete ===\n");
    return 0;
}
