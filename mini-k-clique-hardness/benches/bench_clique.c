/**
 * bench_clique.c ? Performance benchmark for k-Clique algorithms
 *
 * Benchmarks brute-force, color-coding, and Bron-Kerbosch
 * on random graphs of varying sizes.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "../include/kclique_types.h"
#include "../include/kclique_core.h"
#include "../include/kclique_algorithm.h"
#include "../include/kclique_graph.h"

static double get_time_ms(void) {
    return (double)clock() / (CLOCKS_PER_SEC / 1000.0);
}

int main(void) {
    printf("=== k-Clique Algorithm Benchmarks ===\n\n");

    int graph_sizes[] = {10, 15, 20, 30};
    int num_sizes = sizeof(graph_sizes) / sizeof(graph_sizes[0]);

    printf("Brute-force k-Clique: k=3,4 on Erdos-Renyi G(n, 0.5)\n");
    printf("  n   | k | Found |  BF time(ms) | CC time(ms)\n");
    printf("------+---+-------+--------------+-------------\n");

    for (int si = 0; si < num_sizes; si++) {
        int n = graph_sizes[si];
        er_graph_params_t params = {n, 0.5, 42, false, false, 0};
        graph_t *g = graph_generate_erdos_renyi(&params);
        if (!g) continue;

        for (int k = 3; k <= 4 && k <= n; k++) {
            clique_t result = {0};
            double t0 = get_time_ms();
            bool found_bf = find_k_clique_brute(g, k, &result);
            double t1 = get_time_ms();
            double bf_time = t1 - t0;
            if (found_bf) clique_free(&result);

            result = (clique_t){0};
            t0 = get_time_ms();
            bool found_cc = color_coding_solve(g, k, 20, 42, &result);
            t1 = get_time_ms();
            double cc_time = t1 - t0;
            if (found_cc) clique_free(&result);

            printf("  %3d | %d |  %3s  |  %8.2f   |  %8.2f\n",
                   n, k, found_bf ? "YES" : "NO", bf_time, cc_time);
        }
        graph_destroy(g);
    }

    printf("\nBron-Kerbosch Maximum Clique: random G(n, 0.3)\n");
    printf("  n   | max omega | BK time(ms)\n");
    printf("------+-----------+------------\n");

    int bk_sizes[] = {10, 15, 20, 25, 30};
    for (int si = 0; si < 5; si++) {
        int n = bk_sizes[si];
        er_graph_params_t params = {n, 0.3, 99, false, false, 0};
        graph_t *g = graph_generate_erdos_renyi(&params);
        if (!g) continue;

        clique_t result = {0};
        double t0 = get_time_ms();
        int32_t max_k = max_clique_bron_kerbosch(g, &result);
        double t1 = get_time_ms();

        printf("  %3d |    %2d    |  %8.2f\n", n, max_k, t1 - t0);
        clique_free(&result);
        graph_destroy(g);
    }

    printf("\nDone.\n");
    return 0;
}
