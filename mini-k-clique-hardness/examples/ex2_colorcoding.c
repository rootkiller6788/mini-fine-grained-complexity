/**
 * ex2_colorcoding.c ? Color-Coding Algorithm Demo
 *
 * Demonstrates the Alon-Yuster-Zwick color-coding FPT algorithm
 * for k-Clique, showing success probability vs iterations.
 *
 * Knowledge: L5 Algorithms ? Color-coding FPT algorithm
 */

#include <stdio.h>
#include <stdlib.h>
#include "../include/kclique_types.h"
#include "../include/kclique_core.h"
#include <math.h>
#include "../include/kclique_algorithm.h"

int main(void) {
    printf("=== Color-Coding Algorithm Demo ===\n\n");
    printf("Alon, Yuster, Zwick (JACM 1995)\n");
    printf("FPT algorithm: O(2^k * n^2) per coloring,\n");
    printf("success probability >= k!/k^k >= e^{-k} per trial.\n\n");

    /* Create a graph with a hidden 4-clique in 12 vertices */
    int32_t n = 12, k = 4;
    graph_t *g = graph_create(n);

    /* Plant a 4-clique on vertices {2, 5, 8, 11} */
    graph_add_edge(g, 2, 5); graph_add_edge(g, 2, 8); graph_add_edge(g, 2, 11);
    graph_add_edge(g, 5, 8); graph_add_edge(g, 5, 11);
    graph_add_edge(g, 8, 11);

    /* Add random-looking edges that don't create larger cliques */
    graph_add_edge(g, 0, 1); graph_add_edge(g, 0, 3); graph_add_edge(g, 0, 4);
    graph_add_edge(g, 1, 3); graph_add_edge(g, 1, 6); graph_add_edge(g, 1, 7);
    graph_add_edge(g, 3, 4); graph_add_edge(g, 3, 6); graph_add_edge(g, 3, 9);
    graph_add_edge(g, 4, 7); graph_add_edge(g, 4, 9); graph_add_edge(g, 4, 10);
    graph_add_edge(g, 6, 7); graph_add_edge(g, 6, 10);
    graph_add_edge(g, 7, 9); graph_add_edge(g, 7, 10);
    graph_add_edge(g, 9, 10);
    /* Note: none of {0,1,3,4,6,7,9,10} are adjacent to the planted clique
       vertices {2,5,8,11} except within the clique */

    printf("Graph: %d vertices with a planted %d-clique on {2,5,8,11}\n\n", n, k);

    /* Theoretical success probability per trial */
    double p_success = 1.0;
    for (int i = 1; i <= k; i++) p_success *= (double)i / (double)k;
    printf("Theoretical Pr[success per trial] = k!/k^k = %.4f\n", p_success);
    printf("Expected trials needed: 1/p = %.1f\n", 1.0 / p_success);
    printf("e^k = %.1f\n\n", exp((double)k));

    /* Test with increasing iteration counts */
    int iterations[] = {1, 5, 10, 20, 50, 100};
    int num_tests = sizeof(iterations) / sizeof(iterations[0]);

    printf("Iterations | Found? | Clique Vertices\n");
    printf("-----------+--------+----------------\n");
    for (int t = 0; t < num_tests; t++) {
        clique_t result = {0};
        bool found = color_coding_solve(g, k, iterations[t], 42 + t, &result);
        printf("  %8d |  %3s   | ", iterations[t], found ? "YES" : "NO");
        if (found) {
            printf("{");
            for (int i = 0; i < result.k; i++) {
                printf("%d%s", result.vertices[i],
                       i < result.k - 1 ? ", " : "");
            }
            printf("}");
            clique_free(&result);
        } else {
            printf("(not found)");
        }
        printf("\n");
    }

    /* Verify the planted clique is correct */
    printf("\nVerification: ");
    int32_t planted[4] = {2, 5, 8, 11};
    if (is_k_clique(g, planted, 4)) {
        printf("Planting verified (2,5,8,11 form a 4-clique)\n");
    } else {
        printf("ERROR: planting failed!\n");
    }

    graph_destroy(g);
    printf("\nDone.\n");
    return 0;
}
