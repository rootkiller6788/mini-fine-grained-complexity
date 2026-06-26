/**
 * kclique_algorithm.h ? Algorithm declarations for k-Clique
 *
 * FPT algorithms, branch-and-bound, color-coding, and
 * matrix-multiplication approaches for the k-Clique problem.
 *
 * Reference: Alon, Yuster, Zwick, "Color-Coding" (JACM 1995)
 * Reference: Nesetril & Poljak, "Complexity of the subgraph problem" (1985)
 * Reference: Eisenbrand & Grandoni, "Clique in FPT" (2005)
 * Reference: Vassilevska Williams et al., "K-clique via Matrix Mult" (2015)
 *
 * Knowledge: L5 Algorithms ? Color-coding, B&B, matrix approach
 * Knowledge: L8 Advanced Topics ? Derandomization, fast matrix multiplication
 */

#ifndef KCLIQUE_ALGORITHM_H
#define KCLIQUE_ALGORITHM_H

#include "kclique_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ================================================================
 * Color-Coding (Alon, Yuster, Zwick 1995)
 *
 * Randomized FPT algorithm for k-Clique:
 * 1. Randomly color vertices with k colors
 * 2. Use DP to find a k-clique where all vertices have distinct colors
 * 3. If a k-clique exists, Pr[success] >= k!/k^k >= e^{-k}
 * 4. Repeat e^k times for constant success probability
 *
 * Runtime: O((2e)^k * n^2) with randomization
 * Derandomized: O(2^{O(k)} * n^2 log n) using k-perfect hash families
 * ================================================================ */

/**
 * Initialize a color-coding solver state.
 *
 * @param n     number of vertices in the graph
 * @param k     target clique size
 * @param seed  random seed (0 for time-based)
 * @return      initialized color_coding_t, or NULL on failure
 */
color_coding_t *color_coding_init(int32_t n, int32_t k, int64_t seed);

/**
 * Free a color-coding solver state.
 */
void color_coding_destroy(color_coding_t *cc);

/**
 * Assign random colors to all vertices.
 * Each vertex gets a color in {0, ..., k-1} uniformly at random.
 */
void color_coding_random_colors(color_coding_t *cc);

/**
 * Find a k-clique using the color-coding DP approach.
 *
 * DP table: dp[S][v] = true if there exists a clique in vertices
 * colored with the color set S (where |S| = i), ending at vertex v,
 * and all vertices in the clique have distinct colors.
 *
 * Complexity: O(2^k * n^2) time, O(2^k * n) space.
 *
 * @param g     input graph
 * @param cc    color-coding state (with colors already assigned)
 * @param[out] result  found k-clique if successful
 * @return      true if a k-clique was found in this coloring
 */
bool color_coding_find_clique(const graph_t *g,
                               color_coding_t *cc,
                               clique_t *result);

/**
 * Full color-coding k-Clique solver with repetitions.
 * Repeats the random coloring process num_iterations times
 * to achieve high success probability.
 *
 * Success probability after T iterations: 1 - (1 - k!/k^k)^T.
 * With T = e^k, Pr[success] >= 1 - 1/e ~ 0.632.
 * With T = c * e^k, Pr[success] >= 1 - e^{-c}.
 *
 * @param g              input graph
 * @param k              target clique size
 * @param num_iterations number of coloring iterations
 * @param seed           random seed
 * @param[out] result    found clique (if any)
 * @return               true if a k-clique was found
 */
bool color_coding_solve(const graph_t *g, int32_t k,
                         int32_t num_iterations, int64_t seed,
                         clique_t *result);

/**
 * Derandomized color-coding using k-perfect hash families.
 * A k-perfect hash family H of functions [n] -> [k] has the property
 * that for any subset S of [n] of size k, there exists h in H
 * that is injective on S (i.e., assigns distinct colors to S).
 *
 * Using the explicit construction of Naor-Schulman-Srinivasan (1995)
 * or the improved Alon-Gutner (2008) scheme.
 *
 * @param g     input graph
 * @param k     target clique size
 * @param[out] result  found clique
 * @return      true if a k-clique was found
 */
bool color_coding_derandomized(const graph_t *g, int32_t k,
                                clique_t *result);

/* ================================================================
 * Branch-and-Bound for Maximum Clique
 *
 * Carraghan and Pardalos (1990): branch and bound for MCP.
 * Tomita et al. (2003, 2010, 2016): improved B&B with coloring bound.
 *
 * Key idea:
 * - Maintain current best clique
 * - Branch: add one vertex at a time, maintaining clique property
 * - Bound: use greedy graph coloring to compute an upper bound
 *          on the size of a clique in the remaining subgraph
 * - Prune if current + bound <= best
 * ================================================================ */

/**
 * Initialize branch-and-bound solver.
 *
 * @param n         number of vertices
 * @param max_k     maximum clique size to search for
 * @return          initialized state, or NULL
 */
branch_and_bound_t *branch_bound_init(int32_t n, int32_t max_k);

/**
 * Free branch-and-bound state.
 */
void branch_bound_destroy(branch_and_bound_t *bb);

/**
 * Run branch-and-bound to find a k-clique (or maximum clique).
 *
 * Implements depth-first search with greedy-coloring bound.
 * For each recursive step, select candidate vertices that are
 * adjacent to all vertices in the current solution, then apply
 * coloring bound to the candidate set.
 *
 * Complexity: worst-case O(1.442^n) for B&B alone,
 * but much faster in practice due to pruning.
 *
 * @param g     input graph
 * @param k     target clique size (or 0 for maximum)
 * @param bb    initialized branch-and-bound state
 * @param[out] result  found clique
 * @return      size of clique found, or 0 if none >= k
 */
int32_t branch_bound_solve(const graph_t *g, int32_t k,
                            branch_and_bound_t *bb,
                            clique_t *result);

/**
 * Compute a greedy graph coloring of a vertex set.
 * Used as an upper bound: if a graph can be colored with c colors,
 * then any clique has size at most c.
 *
 * Implements the DSATUR (Degree of Saturation) heuristic.
 *
 * @param g              graph
 * @param vertices       subset of vertices to color
 * @param num_vertices   size of the vertex subset
 * @param[out] colors    assigned colors (caller must allocate, size num_vertices)
 * @return               number of colors used
 */
int32_t greedy_coloring_bound(const graph_t *g,
                               const int32_t *vertices,
                               int32_t num_vertices,
                               int32_t *colors);

/**
 * Check if the greedy coloring bound can prune the search.
 *
 * @param current_size  size of current partial clique
 * @param g             graph
 * @param candidates    candidate vertices
 * @param num_candidates  number of candidates
 * @param best_so_far   best clique size found so far
 * @return              true if branch can be pruned
 */
bool coloring_bound_prune(int32_t current_size,
                           const graph_t *g,
                           const int32_t *candidates,
                           int32_t num_candidates,
                           int32_t best_so_far);

/* ================================================================
 * Matrix Multiplication Approach
 *
 * Nesetril & Poljak (1985): k-Clique reducible to triangle detection
 * in a larger graph. For k = 3t, k-Clique reduces to triangle
 * detection in an O(n^{t}) vertex graph.
 *
 * With fast matrix multiplication (omega < 2.373):
 *   - 3-Clique (triangle): O(n^omega) <= O(n^{2.373})
 *   - 4-Clique: O(n^{omega+1}) <= O(n^{3.373})
 *   - 6-Clique: O(n^{3 omega}) <= O(n^{7.119})
 *
 * General: k-Clique in O(n^{omega * ceil(k/3) + k mod 3}).
 * This beats O(n^k) for small values of k.
 * ================================================================ */

/**
 * Detect a triangle (3-clique) using matrix multiplication.
 * Compute A^2 and check if A[i][j] * (A^2)[i][j] > 0 for some i,j.
 * If so, vertices {i, j, k} form a triangle where k is found
 * via the witness of (A^2)[i][j].
 *
 * Complexity: O(n^omega) where omega is the matrix multiplication exponent.
 *
 * @param g        graph (adjacency matrix representation)
 * @param[out] result  triangle vertices if found
 * @return         true if a triangle exists
 */
bool matrix_triangle_detect(const graph_t *g, clique_t *result);

/**
 * Find a 4-clique using matrix multiplication.
 * Partition vertices, check triangles in bipartite graphs,
 * combine to find 4-clique.
 *
 * Complexity: O(n^{omega+1}) = O(n^{3.373}) with current omega.
 *
 * @param g        graph
 * @param[out] result  4-clique vertices
 * @return         true if a 4-clique exists
 */
bool matrix_4clique_detect(const graph_t *g, clique_t *result);

/**
 * General k-clique via matrix multiplication (Nesetril-Poljak method).
 * Reduces k-clique to triangle detection in an auxiliary graph
 * with O(n^{ceil(k/3)}) vertices.
 *
 * Complexity: O(n^{omega * ceil(k/3) + k mod 3}).
 *
 * @param g        graph
 * @param k        clique size
 * @param[out] result  k-clique vertices
 * @return         true if a k-clique exists
 */
bool matrix_kclique_detect(const graph_t *g, int32_t k, clique_t *result);

/**
 * Compute the matrix multiplication exponent omega.
 * Returns the currently best known upper bound (approx 2.371866).
 * Used to estimate algorithm runtimes.
 *
 * @return  best known omega bound
 */
double current_omega_bound(void);

#ifdef __cplusplus
}
#endif

#endif /* KCLIQUE_ALGORITHM_H */
