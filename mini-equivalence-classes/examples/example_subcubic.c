/* Example: Subcubic Equivalence -- APSP, Radius, Diameter, Negative Triangle
 *
 * Demonstrates the subcubic equivalence class by computing APSP via Floyd-Warshall,
 * then extracting all graph metrics (radius, diameter, median) and detecting
 * negative triangles.
 *
 * This example shows:
 *   1. APSP computation on a weighted directed graph
 *   2. Extraction of radius, diameter, median from APSP distances
 *   3. Negative triangle detection via min-plus product
 *   4. Subcubic equivalence verification
 *
 * References:
 *   Floyd (1962), Warshall (1962)
 *   Williams & Williams (2013), STOC
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "subcubic.h"

int main(void) {
    printf("=== Example: Subcubic Equivalence Class ===\n\n");

    /* Build a simple 5-vertex directed graph */
    int32_t n = 5;
    apsp_instance_t *apsp = apsp_create(n);

    /* Define edges: a small directed graph */
    apsp_set_edge(apsp, 0, 1, 3.0);
    apsp_set_edge(apsp, 0, 2, 8.0);
    apsp_set_edge(apsp, 1, 2, 2.0);
    apsp_set_edge(apsp, 1, 3, 5.0);
    apsp_set_edge(apsp, 2, 3, 1.0);
    apsp_set_edge(apsp, 3, 4, 4.0);
    apsp_set_edge(apsp, 2, 4, 7.0);
    apsp_set_edge(apsp, 0, 4, 10.0);

    printf("Graph: 5 vertices, 8 directed edges\n");
    printf("Edges: 0->1(3), 0->2(8), 1->2(2), 1->3(5), 2->3(1), 3->4(4), 2->4(7), 0->4(10)\n\n");

    /* Compute APSP via Floyd-Warshall */
    printf("--- APSP via Floyd-Warshall ---\n");
    bool ok = apsp_floyd_warshall(apsp);
    if (!ok) { printf("Error: negative cycle detected!\n"); apsp_free(apsp); return 1; }

    printf("Shortest path distances:\n");
    for (int32_t i = 0; i < n; i++) {
        for (int32_t j = 0; j < n; j++) {
            double d = apsp_get_distance(apsp, i, j);
            if (d > 1e10) printf("%5s ", "INF");
            else printf("%5.1f ", d);
        }
        printf("\n");
    }
    printf("\n");

    /* Extract graph metrics */
    printf("--- Graph Metrics ---\n");
    double *dist = (double *)malloc((size_t)(n * n) * sizeof(double));
    for (int32_t i = 0; i < n; i++)
        for (int32_t j = 0; j < n; j++)
            dist[i*n+j] = apsp_get_distance(apsp, i, j);

    graph_metrics_t *gm = graph_metrics_create(n, dist);
    graph_metrics_compute(gm);

    printf("Radius:    %.1f (center: vertex %d)\n", gm->radius, gm->center);
    printf("Diameter:  %.1f (path %d -> %d)\n", gm->diameter, gm->diam_u, gm->diam_v);
    printf("Median:    %.1f (median vertex: %d)\n", gm->median, gm->median_vertex);
    printf("Eccentricities: ");
    for (int32_t v = 0; v < n; v++) printf("%.1f ", gm->eccentricity[v]);
    printf("\nUnique center: %s\n\n", graph_metrics_unique_center(gm) ? "yes" : "no");

    /* Negative Triangle detection via Min-Plus product */
    printf("--- Negative Triangle Detection ---\n");
    neg_triangle_instance_t *nt = neg_triangle_create(n);
    neg_triangle_set_edge(nt, 0, 1, 3.0);
    neg_triangle_set_edge(nt, 1, 2, 2.0);
    neg_triangle_set_edge(nt, 2, 0, -6.0);  /* This creates a negative triangle 0->1->2->0 = -1 */

    bool has_neg = neg_triangle_via_min_plus(nt);
    printf("Negative triangle: %s\n", has_neg ? "YES (found)" : "NO (not found)");
    if (has_neg) printf("  Vertices: %d -> %d -> %d\n", nt->vertex_u, nt->vertex_v, nt->vertex_w);

    /* Verify subcubic completeness */
    printf("\n--- Subcubic Completeness ---\n");
    subcubic_status_report();

    /* Cleanup */
    graph_metrics_free(gm);
    free(dist);
    neg_triangle_free(nt);
    apsp_free(apsp);

    printf("\n=== Example Complete ===\n");
    return 0;
}