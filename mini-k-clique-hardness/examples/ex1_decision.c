/**
 * ex1_decision.c ? Decision problem example: k-Clique detection
 *
 * Demonstrates the k-Clique decision problem on a small graph
 * using brute-force and color-coding algorithms.
 *
 * Examples:
 * - Graph with a known 3-clique (triangle)
 * - Graph without a k-clique
 * - Comparison of brute-force vs color-coding
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../include/kclique_types.h"
#include "../include/kclique_core.h"
#include "../include/kclique_algorithm.h"

int main(void) {
    printf("=== k-Clique Decision Problem Example ===\n\n");

    /* Build a graph with a known triangle */
    int n = 8;
    graph_t *g = graph_create(n);
    if (!g) { fprintf(stderr, "Failed to create graph\n"); return 1; }

    /* Triangle on vertices {0, 1, 2} */
    graph_add_edge(g, 0, 1);
    graph_add_edge(g, 0, 2);
    graph_add_edge(g, 1, 2);

    /* Edge (3,4) and star from 3 */
    graph_add_edge(g, 3, 4);
    graph_add_edge(g, 3, 5);
    graph_add_edge(g, 3, 6);

    /* 4-clique on {4,5,6,7} for testing */
    graph_add_edge(g, 4, 5);
    graph_add_edge(g, 4, 6);
    graph_add_edge(g, 4, 7);
    graph_add_edge(g, 5, 6);
    graph_add_edge(g, 5, 7);
    graph_add_edge(g, 6, 7);

    printf("Graph: %d vertices, edges (manually counted)\n", n);
    printf("  Vertices {0,1,2}: triangle\n");
    printf("  Vertices {4,5,6,7}: 4-clique (K4)\n\n");

    /* Decision problems */
    for (int k = 1; k <= 5; k++) {
        printf("k=%d:\n", k);

        /* Brute force */
        clique_t result_brute = {0};
        bool found_brute = find_k_clique_brute(g, k, &result_brute);
        printf("  Brute force:    %s", found_brute ? "FOUND" : "NOT FOUND");
        if (found_brute) {
            printf(" (vertices: {");
            for (int i = 0; i < result_brute.k; i++) {
                printf("%d%s", result_brute.vertices[i],
                       i < result_brute.k - 1 ? ", " : "");
            }
            printf("})");
            clique_free(&result_brute);
        }
        printf("\n");

        /* Color-coding (fewer iterations for demo) */
        clique_t result_cc = {0};
        bool found_cc = color_coding_solve(g, k, 50, 12345, &result_cc);
        printf("  Color-coding:   %s", found_cc ? "FOUND" : "NOT FOUND");
        if (found_cc) {
            printf(" (vertices: {");
            for (int i = 0; i < result_cc.k; i++) {
                printf("%d%s", result_cc.vertices[i],
                       i < result_cc.k - 1 ? ", " : "");
            }
            printf("})");
            clique_free(&result_cc);
        }
        printf("\n\n");
    }

    /* Count cliques */
    printf("Clique counts in this graph:\n");
    for (int k = 1; k <= 4; k++) {
        int64_t cnt = count_k_cliques(g, k);
        printf("  %d-cliques: %lld\n", k, (long long)cnt);
    }

    graph_destroy(g);
    printf("\nDone.\n");
    return 0;
}
