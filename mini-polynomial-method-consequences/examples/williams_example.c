#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "poly_method.h"

/* ============================================================================
 * williams_example.c - Williams' Algorithm via Polynomial Method
 *
 * L5: Algorithm - Williams-style SAT solving via polynomial method
 * L6: Canonical problem - CNF-SAT
 * L7: Application - circuit lower bounds connection
 * ============================================================================ */

int main(void) {
    printf("=== Williams' Algorithm: Polynomial Method for SAT ===\n\n");

    /* Create a 3-SAT formula with n=8 variables, m=6 clauses */
    int32_t n = 8, m = 6, k = 3;
    printf("Creating random 3-SAT instance: n=%d, m=%d, k=%d\n", n, m, k);

    williams_cnf_t *f = williams_generate_ksat(n, m, k, 42);
    if (!f) { printf("Failed to create formula\n"); return 1; }

    printf("Clauses:\n");
    int32_t clause_size = k + 1;
    for (int32_t c = 0; c < f->num_clauses; c++) {
        printf("  C%d: (", c + 1);
        int32_t base = c * clause_size;
        for (int32_t i = 0; f->clause_data[base + i] != 0; i++) {
            int32_t lit = f->clause_data[base + i];
            if (i > 0) printf(" OR ");
            if (lit < 0) printf("~");
            printf("x%d", (lit > 0) ? lit : -lit);
        }
        printf(")\n");
    }

    /* Brute force baseline */
    printf("\n--- Brute Force (O(m * 2^n)) ---\n");
    williams_result_t r_bf = williams_brute_force(f);
    printf("SAT: %s\n", r_bf.satisfiable ? "YES" : "NO");
    printf("Time: %.3f ms\n", r_bf.runtime_ms);
    printf("Evaluations: %llu\n", (unsigned long long)r_bf.evaluations);
    if (r_bf.satisfiable) {
        printf("Witness: [");
        for (int32_t i = 0; i < n; i++) printf("%d ", r_bf.witness[i]);
        printf("]\n");
        printf("Verified: %s\n", williams_verify_witness(f, r_bf.witness) ? "PASS" : "FAIL");
    }

    /* Williams polynomial method */
    printf("\n--- Williams Polynomial Method (O*(2^{n/2})) ---\n");
    williams_result_t r_w = williams_solve(f);
    printf("SAT: %s\n", r_w.satisfiable ? "YES" : "NO");
    printf("Time: %.3f ms\n", r_w.runtime_ms);
    printf("Evaluations: %llu\n", (unsigned long long)r_w.evaluations);
    if (r_w.satisfiable) {
        printf("Witness: [");
        for (int32_t i = 0; i < n; i++) printf("%d ", r_w.witness[i]);
        printf("]\n");
        printf("Verified: %s\n", williams_verify_witness(f, r_w.witness) ? "PASS" : "FAIL");
    }

    /* Demonstrate residual polynomial construction */
    printf("\n--- Residual Polynomial Demo ---\n");
    int32_t split = n / 2;
    williams_assign_t *partial = williams_assign_create(n);
    /* Fix first 4 variables to a specific assignment */
    for (int32_t i = 0; i < split; i++)
        partial->values[i] = (int8_t)(i % 2);  /* alternating */

    printf("Partial assignment (vars 1-%d): [", split);
    for (int32_t i = 0; i < split; i++)
        printf("%d ", partial->values[i]);
    printf("]\n");

    polynomial_t *residual = williams_residual_polynomial(f, partial, split);
    if (residual) {
        printf("Residual polynomial: vars=%d, terms=%d, degree=%d\n",
               residual->num_vars, residual->num_terms, residual->degree);

        /* Show some terms */
        int32_t to_show = (residual->num_terms < 10) ? residual->num_terms : 10;
        printf("First %d terms of residual polynomial:\n", to_show);
        for (int32_t i = 0; i < to_show; i++) {
            printf("  Term %d: coeff=%d, vars=0x%llx, deg=%d\n",
                   i, (int)residual->terms[i].coeff & 1,
                   (unsigned long long)residual->terms[i].mono.vars,
                   residual->terms[i].mono.degree);
        }
        if (residual->num_terms > 10)
            printf("  ... (%d more terms)\n", residual->num_terms - 10);

        poly_destroy(residual);
    }
    williams_assign_destroy(partial);

    /* Comparison summary */
    printf("\n--- Comparison ---\n");
    printf("Brute force time: %.3f ms\n", r_bf.runtime_ms);
    printf("Williams time:    %.3f ms\n", r_w.runtime_ms);
    if (r_bf.satisfiable == r_w.satisfiable)
        printf("Results match: YES\n");
    else
        printf("Results match: NO (MISMATCH!)\n");

    williams_cnf_destroy(f);
    free(r_bf.witness);
    free(r_w.witness);

    printf("\n=== Example complete ===\n");
    return 0;
}
