#ifndef APSP_H
#define APSP_H
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdbool.h>
#include <float.h>
#ifdef __cplusplus
extern "C" {
#endif
#define APSP_INF DBL_MAX
typedef struct { double *dist; int32_t *next; int32_t n; bool has_negative_cycle; } apsp_graph_t;
typedef struct { double *weights; int32_t n; } apsp_adjacency_t;
typedef struct { int32_t *vertices; int32_t length; double total_weight; } apsp_path_t;
typedef enum { APSP_HARD_BASELINE=0,APSP_HARD_NEG_TRIANGLE=1,APSP_HARD_MIN_PLUS_PRODUCT=2,APSP_HARD_RADIUS=3,APSP_HARD_MEDIAN=4,APSP_HARD_BETWEENNESS=5,APSP_HARD_COUNT=6 } apsp_hardness_class_t;
typedef struct { const char *statement; double epsilon_min; const char *history; } apsp_conjecture_t;
typedef struct { double *data; int32_t rows; int32_t cols; } apsp_matrix_t;
apsp_graph_t *apsp_floyd_warshall(const apsp_adjacency_t *adj);
apsp_graph_t *apsp_johnson(apsp_adjacency_t *adj);
apsp_graph_t *apsp_dijkstra_all(const apsp_adjacency_t *adj);
bool apsp_has_negative_cycle(const apsp_adjacency_t *adj);
apsp_path_t *apsp_extract_path(const apsp_graph_t *g, int32_t u, int32_t v);
double apsp_graph_diameter(const apsp_graph_t *g);
apsp_matrix_t *apsp_min_plus_multiply(const apsp_matrix_t *A, const apsp_matrix_t *B);
apsp_conjecture_t apsp_conjecture_status(void);
apsp_adjacency_t *apsp_adjacency_create(int32_t n);
void apsp_adjacency_destroy(apsp_adjacency_t *adj);
void apsp_adjacency_set_edge(apsp_adjacency_t *adj, int32_t u, int32_t v, double w);
apsp_graph_t *apsp_graph_create(int32_t n);
void apsp_graph_destroy(apsp_graph_t *g);
void apsp_print_graph(const apsp_graph_t *g);
apsp_adjacency_t *apsp_generate_random_graph(int32_t n, double density, double max_weight, uint64_t seed);
#ifdef __cplusplus
}
#endif
#endif

/*
 * ============================================================
 * APSP: Extended Documentation
 * ============================================================
 *
 * The All-Pairs Shortest Paths problem is one of the most
 * fundamental problems in graph algorithms. Given a weighted
 * directed graph with n vertices and m edges, compute the
 * shortest-path distance between every pair of vertices.
 *
 * Classic algorithms:
 *   Floyd-Warshall (1962): O(n^3) time, O(n^2) space.
 *     Dynamic programming over intermediate vertices.
 *   Johnson (1977): O(nm + n^2 log n) for sparse graphs.
 *     Uses Bellman-Ford for reweighting + n x Dijkstra.
 *   Dijkstra (1959): O(m + n log n) per source with Fibonacci
 *     heaps, giving O(nm + n^2 log n) overall.
 *
 * The APSP Conjecture (Vassilevska-Williams 2009):
 *   For every epsilon > 0, there is no algorithm that solves
 *   APSP on n-vertex graphs with arbitrary real edge weights
 *   in O(n^{3-epsilon}) time.
 *
 * Subcubic equivalences:
 *   APSP is subcubic-equivalent to:
 *   - Min-Plus Matrix Product
 *   - Negative Triangle Detection
 *   - Graph Radius
 *   - Graph Median
 *   - Betweenness Centrality
 *   - Diameter (in weighted graphs)
 *
 *   This means: if ANY of these problems can be solved in
 *   O(n^{3-epsilon}) time, then ALL of them can.
 *
 * Breakthrough results:
 *   - Seidel (1995): APSP for unweighted undirected graphs
 *     in O(M(n) log n) using fast matrix multiplication,
 *     where M(n) = O(n^{2.3729}) is the matrix mult exponent.
 *   - Zwick (2002): O(n^{2.575}) for directed graphs with
 *     small integer weights.
 *   - Williams (2014): O(n^3 / 2^{Omega(sqrt(log n))}) for
 *     general real weights, the first algorithm to beat n^3
 *     by more than a polylogarithmic factor.
 *
 * Despite these results, the O(n^{3-epsilon}) barrier remains:
 *   no known algorithm achieves O(n^{2.999}) for general
 *   real-weighted graphs.
 *
 * The Min-Plus Semiring:
 *   APSP is intimately connected to the min-plus (tropical)
 *   semiring (R U {inf}, min, +). In this semiring:
 *   - Matrix "multiplication" is: C[i][j] = min_k{A[i][k]+B[k][j]}
 *   - APSP = (n-1)-fold min-plus power of adjacency matrix
 *   - Repeated squaring: D = A^{(n-1)} = A^{(2^{ceil(log n)})}
 *     requires O(log n) min-plus products
 *
 *   However, unlike standard matrix multiplication over a ring,
 *   min-plus product currently has no known O(n^{3-epsilon})
 *   algorithm, despite being the direct analog.
 *
 * For further reading:
 *   - Floyd (1962): "Algorithm 97: Shortest Path", CACM
 *   - Warshall (1962): "A theorem on Boolean matrices", JACM
 *   - Johnson (1977): "Efficient algorithms for shortest paths
 *     in sparse networks", JACM
 *   - Vassilevska-Williams & Williams (2010): "Subcubic
 *     Equivalences...", FOCS 2010
 */
/*
 * APSP Algorithm Comparison:
 *
 * | Algorithm          | Time               | Space  | Neg Edges |
 * |--------------------|--------------------|--------|-----------|
 * | Floyd-Warshall     | Theta(n^3)         | O(n^2) | Yes       |
 * | Johnson (1977)     | O(nm + n^2 log n)  | O(n^2) | Yes       |
 * | Dijkstra-All       | O(nm + n^2 log n)  | O(n^2) | No        |
 * | Bellman-Ford-All   | O(n^2 m)           | O(n^2) | Yes       |
 * | Rep. Squaring      | O(M(n) log n)      | O(n^2) | Yes       |
 * | Williams (2014)    | O(n^3/2^{sqrt(log n)}) | O(n^2) | Yes  |
 *
 * where M(n) = time for one min-plus matrix product.
 *
 * For dense graphs (m = Theta(n^2)):
 *   - Floyd-Warshall:  O(n^3)
 *   - Johnson:         O(n^3) (Dijkstra dominates)
 *   - Williams 2014:   O(n^3 / 2^{sqrt(log n)})  [best known]
 *
 * For sparse graphs (m = O(n)):
 *   - Johnson:         O(n^2 log n)
 *   - Dijkstra-All:    O(n^2 log n) [no negative edges]
 */