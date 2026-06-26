#ifndef SUBCUBIC_H
#define SUBCUBIC_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <math.h>
#include <float.h>

/*
 * subcubic.h -- Subcubic Equivalence Class
 *
 * The subcubic equivalence class contains problems whose best known
 * algorithms run in O(n^3) time and for which no O(n^{3-epsilon})
 * algorithm is known for any epsilon > 0.
 *
 * Theorem (Williams & Williams 2013): Under subcubic reductions, the
 * following problems are equivalent:
 *   - All-Pairs Shortest Paths (APSP)
 *   - Negative Triangle detection
 *   - Min-Plus Matrix Multiplication
 *   - Radius, Median, Diameter
 *   - Second Shortest Path
 *   - Replacement Paths
 *
 * Canonical problem: APSP (or equivalently, Min-Plus product)
 * Conjecture: No O(n^{3-epsilon}) for APSP in dense graphs.
 *
 * Best known: O(n^3/log^2 n) (Williams 2014) or O(n^omega)
 * with omega < 2.372 (Alman & Williams 2020).
 *
 * L1: APSP, Negative Triangle, Radius/Median/Diameter definitions
 * L2: Subcubic reduction, subcubic equivalence concept
 * L3: Graph adjacency matrices, distance matrices (min-plus semiring)
 * L4: Williams & Williams (2013) subcubic equivalence theorem
 * L5: Floyd-Warshall, Min-Plus product, fast matrix multiplication
 * L6: APSP, Negative Triangle, Radius, Median, Diameter
 */

/* ---- L1: APSP Instance ---- */

typedef struct {
    int32_t     n;
    double     *dist;
    double     *next_dist;
    bool        has_negative_cycle;
    int32_t     negative_cycle_vertex;
} apsp_instance_t;

/* ---- L1: Negative Triangle Instance ---- */

typedef struct {
    int32_t     n;
    double     *weight;
    bool        found;
    int32_t     vertex_u;
    int32_t     vertex_v;
    int32_t     vertex_w;
} neg_triangle_instance_t;

/* ---- L1: Graph Metrics (Radius, Median, Diameter) ---- */

typedef struct {
    int32_t     n;
    double     *dist;
    double      radius;
    int32_t     center;
    double      diameter;
    int32_t     diam_u;
    int32_t     diam_v;
    double      median;
    int32_t     median_vertex;
    double     *eccentricity;
} graph_metrics_t;

/* ---- L1: Min-Plus Matrix Product ---- */

typedef struct {
    int32_t     n;
    double     *A;
    double     *B;
    double     *C;
} min_plus_product_t;

/* ---- L2/L5: APSP Operations ---- */

apsp_instance_t *apsp_create(int32_t n);
void apsp_free(apsp_instance_t *inst);
void apsp_set_edge(apsp_instance_t *inst, int32_t u, int32_t v, double weight);
double apsp_get_edge(const apsp_instance_t *inst, int32_t u, int32_t v);
double apsp_get_distance(const apsp_instance_t *inst, int32_t u, int32_t v);

/*
 * Floyd-Warshall algorithm (1962).
 * Classic DP for APSP. Time: Theta(n^3), Space: Theta(n^2).
 * Handles negative edge weights but not negative cycles.
 */
bool apsp_floyd_warshall(apsp_instance_t *inst);

/*
 * APSP via repeated squaring of the min-plus matrix product.
 * Establishes APSP <=_subcubic Min-Plus Product.
 * Time: O(T(n) log n) where T(n) = min-plus multiplication time.
 */
bool apsp_via_min_plus(apsp_instance_t *inst, bool (*min_plus)(min_plus_product_t *));

/*
 * Floyd-Warshall path reconstruction.
 * Returns number of vertices in the path, 0 if no path, -1 if negative cycle.
 */
int32_t apsp_reconstruct_path(const apsp_instance_t *inst, int32_t u, int32_t v,
                               int32_t *path, int32_t max_path_len);

/* ---- L5: Negative Triangle Detection ---- */

neg_triangle_instance_t *neg_triangle_create(int32_t n);
void neg_triangle_free(neg_triangle_instance_t *inst);
void neg_triangle_set_edge(neg_triangle_instance_t *inst, int32_t u, int32_t v, double w);

/*
 * Negative Triangle via APSP.
 * Theorem: NegTriangle <=_subcubic APSP.
 * Run APSP and check if dist[v][v] < 0 for any v.
 */
bool neg_triangle_via_apsp(neg_triangle_instance_t *inst);

/*
 * Negative Triangle via Min-Plus product.
 * Theorem: NegTriangle <=_subcubic Min-Plus Product.
 * Time: O(n^omega) using fast matrix multiplication.
 */
bool neg_triangle_via_min_plus(neg_triangle_instance_t *inst);

/*
 * APSP via Negative Triangle oracle.
 * Theorem: APSP <=_subcubic NegTriangle.
 * The log-factor blowup shows APSP and NegTriangle are subcubic equivalent.
 */
bool apsp_via_neg_triangle(apsp_instance_t *inst);

/* ---- L5: Graph Metrics ---- */

graph_metrics_t *graph_metrics_create(int32_t n, const double *dist);
void graph_metrics_free(graph_metrics_t *gm);

/*
 * Compute radius, diameter, median, eccentricity from APSP distances.
 * Theorem (Williams & Williams 2013): These are all subcubic-equivalent to APSP.
 * Time: O(n^2) after APSP distances are computed.
 */
void graph_metrics_compute(graph_metrics_t *gm);

/*
 * Check if the graph has a unique center.
 * Returns 1 if exactly one center, 0 if multiple centers.
 */
int graph_metrics_unique_center(const graph_metrics_t *gm);

/* ---- L5: Min-Plus Matrix Product ---- */

min_plus_product_t *min_plus_create(int32_t n);
void min_plus_free(min_plus_product_t *mp);

/*
 * Standard min-plus multiplication: C[i][j] = min_k(A[i][k] + B[k][j]).
 * Time: Theta(n^3) naive.
 */
void min_plus_naive(min_plus_product_t *mp);

/*
 * Fast min-plus multiplication via fast matrix multiplication.
 * Uses (min,+) semiring property. Time: O(n^omega * log M).
 * Reference: Alon, Galil, Margalit (1997), J. Computer and System Sciences.
 */
void min_plus_fast(min_plus_product_t *mp);

/*
 * Verify that two min-plus products are equal (within tolerance).
 */
bool min_plus_verify(const min_plus_product_t *a, const min_plus_product_t *b, double eps);

/* ---- L6: Subcubic Completeness ---- */

/*
 * Verify that a problem P is subcubic-complete:
 * 1. P can be solved in O(n^3) time.
 * 2. APSP reduces to P in subcubic time.
 */
bool verify_subcubic_completeness(int32_t n,
    bool (*solve_P)(void *instance, int32_t n),
    void *(*reduce_APSP_to_P)(const apsp_instance_t *apsp));

/*
 * Print a summary of the subcubic equivalence class.
 */
void subcubic_status_report(void);

#endif /* SUBCUBIC_H */