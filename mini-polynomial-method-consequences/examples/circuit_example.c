#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "poly_method.h"

/* ============================================================================
 * circuit_example.c - Circuit Lower Bounds via Polynomial Method
 *
 * L8: Advanced topic - Razborov-Smolensky lower bounds
 * L9: Research frontier - natural proofs barrier
 * ============================================================================ */

int main(void) {
    printf("=== Circuit Lower Bounds via Polynomial Method ===\n\n");

    /* Demonstrate the polynomial degree gap for MOD functions */
    printf("--- Polynomial Degree of MOD Functions ---\n");
    printf("(Razborov-Smolensky 1987: MOD_q not in AC^0[p] for distinct primes p,q)\n\n");

    for (int32_t n = 2; n <= 8; n += 2) {
        printf("n=%d:\n", n);

        /* MOD_3 over GF(2) */
        int32_t d32 = mod_m_degree_over_p(n, 3, 2);
        printf("  MOD_3 over GF(2): degree = %d", d32);

        /* MOD_2 over GF(2) */
        int32_t d22 = mod_m_degree_over_p(n, 2, 2);
        printf("  |  MOD_2 over GF(2): degree = %d", d22);

        /* MOD_3 over GF(3) */
        int32_t d33 = mod_m_degree_over_p(n, 3, 3);
        printf("  |  MOD_3 over GF(3): degree = %d\n", d33);
    }

    printf("\nInterpretation:\n");
    printf("  MOD_2 over GF(2) has degree 1 (PARITY = x1+x2+...+xn)\n");
    printf("  MOD_3 over GF(3) has degree <= 2 (by Fermat: sum^{p-1})\n");
    printf("  MOD_3 over GF(2) has degree Omega(n) => not in AC^0[2]\n");

    /* Demonstrate MAJORITY degree */
    printf("\n--- MAJORITY Function Degree ---\n");
    for (int32_t n = 2; n <= 8; n += 2) {
        int32_t deg = majority_degree_over_p(n, 2);
        printf("  n=%d: MAJORITY degree over GF(2) = %d (max=%d)\n", n, deg, n);
    }
    printf("  => MAJORITY not in AC^0[p] (Razborov 1987)\n");

    /* Razborov approximation analysis */
    printf("\n--- Razborov Approximation Analysis ---\n");
    razborov_approx_params_t params[] = {
        {2, 100, 4, 2},
        {3, 1000, 8, 2},
        {4, 10000, 12, 2},
    };
    for (int32_t i = 0; i < 3; i++) {
        int32_t bound = razborov_degree_bound(&params[i], 0.1);
        printf("  Depth=%d, Size=%d, Fan-in=%d, p=%d => degree bound = %d\n",
               params[i].depth, params[i].size, params[i].max_fan_in,
               params[i].p, bound);
        printf("    Functions requiring degree > %d cannot be computed\n", bound);
    }

    /* Algebraic circuit example */
    printf("\n--- Algebraic Circuit Example ---\n");
    algebraic_circuit_t *ac = ac_create(20, 4);
    if (ac) {
        /* Build circuit for (x0 + x1) * (x2 + x3) */
        int32_t g0 = ac_add_input(ac, 0);
        int32_t g1 = ac_add_input(ac, 1);
        int32_t g2 = ac_add_input(ac, 2);
        int32_t g3 = ac_add_input(ac, 3);
        int32_t s01 = ac_add_add(ac, g0, g1);
        int32_t s23 = ac_add_add(ac, g2, g3);
        int32_t prod = ac_add_mul(ac, s01, s23);
        ac->output_gate = prod;

        printf("Circuit: (x0 + x1) * (x2 + x3)\n");
        printf("  Gates: %d\n", ac->num_gates);
        printf("  Formal degree: %d\n", ac_compute_degree(ac));

        double inputs[] = {1.0, 0.0, 1.0, 1.0};
        printf("  Input: [%.0f, %.0f, %.0f, %.0f]\n",
               inputs[0], inputs[1], inputs[2], inputs[3]);
        printf("  Output: %.0f (expected: (1+0)*(1+1)=2)\n",
               ac_evaluate(ac, inputs));

        ac_destroy(ac);
    }

    /* Natural proofs barrier discussion */
    printf("\n--- Natural Proofs Barrier ---\n");
    printf("Razborov-Rudich (1997): Most known lower bound techniques\n");
    printf("are 'natural' and cannot prove P != NP under standard\n");
    printf("cryptographic assumptions (existence of PRGs).\n\n");
    printf("The polynomial method:\n");
    printf("  - Uses algebraic structure (polynomials over finite fields)\n");
    printf("  - Leverages specific function properties (MOD_m)\n");
    printf("  - Combines circuit analysis with algorithm design\n");
    printf("  - Can sometimes bypass the natural proofs barrier\n\n");
    printf("Williams (2014): NEXP not in ACC^0, proved via\n");
    printf("polynomial method + algorithm design.\n");

    /* Boolean function analysis */
    printf("\n--- Boolean Function Analysis ---\n");
    int32_t n_small = 3;
    uint8_t *tt_maj = build_majority_truth_table(n_small);
    if (tt_maj) {
        double *inf = compute_influences(tt_maj, n_small);
        if (inf) {
            printf("Influence of variables on MAJORITY(n=3):\n");
            for (int32_t i = 0; i < n_small; i++)
                printf("  Inf[x%d] = %.3f\n", i, inf[i]);
            free(inf);
        }

        bool low_conc = function_has_low_degree_concentration(tt_maj, n_small, 1, 0.2);
        printf("Low-degree concentration (d=1, eps=0.2): %s\n",
               low_conc ? "YES" : "NO");
        free(tt_maj);
    }

    printf("\n=== Example complete ===\n");
    return 0;
}
