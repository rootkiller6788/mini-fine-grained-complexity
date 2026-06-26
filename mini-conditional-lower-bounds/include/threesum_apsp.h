/* threesum_apsp.h -- 3SUM and APSP Conjectures
 *
 * 3SUM Conjecture (Gajentaan-Overmars, 1995; Patrascu, 2010):
 *   In the word-RAM model, 3SUM on n numbers requires n^{2-o(1)} time.
 *
 * The 3SUM problem: Given a set of n integers, determine if there
 * exists a triple a, b, c such that a + b + c = 0.
 *
 * APSP Conjecture (Williams-Williams, 2010):
 *   All-Pairs Shortest Paths on n-vertex graphs requires n^{3-o(1)} time
 *   in the word-RAM model.
 *
 * These two conjectures, together with SETH and OVC, form the four
 * pillars of fine-grained complexity. Each generates its own web
 * of conditional lower bounds.
 *
 * Key 3SUM-hard problems:
 *   3-Points-on-a-Line, Polygon Containment, Minimum Area Triangle,
 *   Motion Planning, etc.
 *
 * Key APSP-hard problems:
 *   Negative Triangle, Radius, Median, Betweenness Centrality,
 *   Replacement Paths, Second Shortest Path, etc.
 *
 * References:
 *   Gajentaan & Overmars (1995): "On a Class of O(n^2) Problems
 *        in Computational Geometry"
 *   Patrascu (2010): "Towards Polynomial Lower Bounds for Dynamic Problems"
 *   Williams & Williams (2010): "Subcubic Equivalences Between
 *        Path, Matrix, and Triangle Problems"
 *   Vassilevska Williams & Williams (2018): "Subcubic Equivalences
 *        Between Path, Matrix, and Triangle Problems" (journal version)
 */

#ifndef THREESUM_APSP_H
#define THREESUM_APSP_H

#include "condlb.h"
#include <stdint.h>
#include <limits.h>

/* ============================================================================
 * L1: Definitions
 * ============================================================================ */

/* 3SUM instance: n integers (can be large, to avoid trivial hashing) */
typedef struct {
    int64_t* numbers;
    int      n;
    int      capacity;
    int      sorted;    /* 1 if numbers are sorted */
} ThreeSumInstance;

/* A 3SUM solution: i, j, k such that A[i] + A[j] + A[k] = 0 */
typedef struct {
    int     i, j, k;    /* indices of the three numbers */
    int64_t sum;        /* should be 0 */
    int     found;      /* 1 if solution exists */
} ThreeSumSolution;

/* APSP problem: weighted directed graph */
typedef struct {
    int      n_vertices;
    int64_t* adj_matrix; /* n*n matrix, INF for no edge */
    int64_t** dist;      /* output: n*n shortest distances */
    int**    next;       /* output: n*n next vertex for path reconstruction */
} ApspInstance;

/* ============================================================================
 * L5: Algorithms - 3SUM Solvers
 * ============================================================================ */

/* Create/destroy 3SUM instance */
ThreeSumInstance* ts3_create(int capacity);
void              ts3_free(ThreeSumInstance* inst);
void              ts3_add(ThreeSumInstance* inst, int64_t x);
void              ts3_sort(ThreeSumInstance* inst);

/* Naive O(n^3) 3SUM: check all triples.
 * Time: O(n^3). */
ThreeSumSolution ts3_brute_force(const ThreeSumInstance* inst);

/* O(n^2) 3SUM using sorting + two-pointer technique.
 * Time: O(n^2). This is the baseline that 3SUM conjecture
 * claims is essentially optimal. */
ThreeSumSolution ts3_quadratic(const ThreeSumInstance* inst);

/* Subquadratic 3SUM using FFT when numbers are small.
 * Time: O(n * log n + U * log U) where U = range of values.
 * Shows that 3SUM is easier for small integers. */
ThreeSumSolution ts3_fft_small_range(const ThreeSumInstance* inst, int max_val);

/* Randomized subquadratic 3SUM using hashing (Baran-Demaine-Patrascu, 2005).
 * Time: O~(n^2 / log^2 n) in expectation.
 * The polylog speedup shows the fine-grained nature of the bound. */
ThreeSumSolution ts3_randomized_hashing(const ThreeSumInstance* inst);

/* 3SUM decision variant: return 1 if any triple sums to 0. */
int ts3_exists(const ThreeSumInstance* inst);

/* Count all 3SUM solutions. Time: O(n^2). */
int ts3_count_all(const ThreeSumInstance* inst);

/* ============================================================================
 * L5: Algorithms - APSP Solvers
 * ============================================================================ */

/* Create/destroy APSP instance */
ApspInstance* apsp_create(int n);
void          apsp_free(ApspInstance* inst);
void          apsp_add_edge(ApspInstance* inst, int u, int v, int64_t w);

/* Floyd-Warshall: O(n^3). The baseline algorithm.
 * Reference: Floyd (1962), Warshall (1962) */
void apsp_floyd_warshall(ApspInstance* inst);

/* Dijkstra from each vertex: O(n * (m + n log n)).
 * Better for sparse graphs. */
void apsp_dijkstra_all(ApspInstance* inst);

/* Repeated squaring / min-plus matrix multiplication approach:
 * O(n^3 log n) with standard MM, O(n^omega) with fast MM.
 * Reference: Fisher-Meyer (1971), Munro (1971) */
void apsp_min_plus_product(ApspInstance* inst);

/* Subcubic APSP for small integer weights via fast MM.
 * When weights are bounded, APSP is in O~(n^(3+omega)/2).
 * Reference: Zwick (2002), Williams (2014) */
void apsp_small_weights(ApspInstance* inst, int max_weight);

/* Check for negative cycles in the graph.
 * Time: O(n^3) using Floyd-Warshall variant.
 * Returns 1 if negative cycle exists. */
int apsp_has_negative_cycle(const ApspInstance* inst);

/* Reconstruct shortest path from u to v.
 * Returns number of vertices on path, -1 if unreachable. */
int apsp_reconstruct_path(const ApspInstance* inst, int u, int v, int* path, int max_len);

/* ============================================================================
 * L4: Theorems - Conditional Lower Bounds
 * ============================================================================ */

/* Theorem (Patrascu, 2010): 3SUM requires Omega(n^2) operations in
 * the 3-linear decision tree model.
 * This function checks if an algorithm's complexity contradicts the
 * 3SUM conjecture. Returns 1 if conjecture would be refuted. */
int ts3_check_conjecture(double algo_exponent);

/* Theorem (Williams-Williams, 2010): APSP and Negative Triangle are
 * subcubic-equivalent: APSP requires n^{3-o(1)} iff
 * Negative Triangle requires n^{3-o(1)}.
 * Compute the subcubic equivalence exponent: if one problem requires
 * n^c, the other also requires n^c. */
double apsp_subcubic_equivalence_exponent(void);

/* Theorem (Vassilevska-Williams-Williams, 2018):
 * The subcubic equivalence class has >30 problems.
 * This function returns the current count of known problems in the
 * subcubic equivalence class. */
int apsp_subcubic_class_size(void);

/* ============================================================================
 * L7: Applications
 * ============================================================================ */

/* 3SUM application: Check if n points in the plane are in general
 * position (no three collinear). Reduces to 3SUM.
 * Time: O(n^2). */
int geom_three_collinear(const int64_t* x, const int64_t* y, int n);

/* APSP application: Compute the graph radius = min_v max_u dist(v,u).
 * Equivalent to APSP under subcubic reductions. */
int64_t apsp_graph_radius(const ApspInstance* inst);

/* APSP application: Betweenness centrality for all vertices.
 * Brute force: O(n^3) using APSP. */
void apsp_betweenness_centrality(const ApspInstance* inst, double* centrality);

/* 3SUM application: Minimum area triangle among n points.
 * Reduces to 3SUM. Time: O(n^2). */
int64_t geom_min_area_triangle(const int64_t* x, const int64_t* y, int n);

#endif /* THREESUM_APSP_H */
