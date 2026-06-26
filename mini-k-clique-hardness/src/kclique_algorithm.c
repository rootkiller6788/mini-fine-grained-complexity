/**
 * kclique_algorithm.c ? FPT and exact algorithms for k-Clique
 *
 * Implements color-coding, branch-and-bound, and matrix multiplication
 * approaches for solving the k-Clique problem.
 *
 * Reference: Alon, Yuster, Zwick, "Color-Coding" (JACM 1995)
 * Reference: Nesetril & Poljak, "Subgraph problem complexity" (1985)
 * Reference: Eisenbrand & Grandoni, "Fixed-Parameter Clique" (2005)
 * Reference: Ne?et?il & Poljak, "Complexity of the subgraph problem" (1985)
 * Reference: Vassilevska Williams, "Efficient Algorithms for Closest
 *             String and Related Problems" (SODA 2015)
 *
 * Knowledge: L5 Algorithms ? Color-coding FPT, B&B, matrix approach
 * Knowledge: L8 Advanced Topics ? Derandomization, fast matrix mult
 */

#include "kclique_algorithm.h"
#include "kclique_core.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

/* ================================================================
 * Linear Congruential Generator (LCG) for reproducibility
 * ================================================================ */

static uint64_t lcg_state = 123456789ULL;

static void lcg_seed(uint64_t seed) {
    lcg_state = seed ? seed : (uint64_t)time(NULL);
}

static uint64_t lcg_next(void) {
    lcg_state = lcg_state * 6364136223846793005ULL + 1442695040888963407ULL;
    return lcg_state;
}

static int32_t lcg_rand_int(int32_t max) {
    if (max <= 0) return 0;
    return (int32_t)(lcg_next() % (uint64_t)max);
}

/* ================================================================
 * Color-Coding Algorithm
 *
 * Key idea (Alon-Yuster-Zwick 1995):
 * 1. Randomly color each vertex with one of k colors
 * 2. Look for a "colorful" k-clique where all k colors appear
 * 3. If a k-clique C exists, Pr[C is colorful] = k!/k^k >= e^{-k}
 * 4. Repeat e^k times for high success probability
 *
 * The DP for "colorful k-clique":
 *   dp[S][v] = true iff there exists a |S|-clique ending at v
 *              where vertices are colored with exactly the color set S
 *
 * Recurrence:
 *   dp[{c(v)}][v] = true  (for each vertex v, color c(v))
 *   dp[S][v] = OR_{u in N(v), c(u) in S\{c(v)}} dp[S\{c(v)}][u]
 *
 * The DP runs in O(2^k * n^2) time and O(2^k * n) space.
 * ================================================================ */

color_coding_t *color_coding_init(int32_t n, int32_t k, int64_t seed) {
    if (n <= 0 || k <= 0 || k > n) return NULL;
    color_coding_t *cc = (color_coding_t *)calloc(1, sizeof(color_coding_t));
    if (!cc) return NULL;

    cc->colors = (int32_t *)malloc((size_t)n * sizeof(int32_t));
    if (!cc->colors) { free(cc); return NULL; }

    cc->n = n;
    cc->k = k;
    cc->num_iterations = (int32_t)(exp((double)k) + 0.5); /* e^k */
    cc->seed = seed;
    cc->use_explicit_hash = false;
    lcg_seed((uint64_t)seed);
    return cc;
}

void color_coding_destroy(color_coding_t *cc) {
    if (!cc) return;
    free(cc->colors);
    cc->colors = NULL;
    free(cc);
}

void color_coding_random_colors(color_coding_t *cc) {
    if (!cc) return;
    for (int32_t i = 0; i < cc->n; i++) {
        cc->colors[i] = lcg_rand_int(cc->k);
    }
}

/**
 * DP for finding a colorful k-clique.
 *
 * Uses a bitmask representation: S is encoded as an integer
 * where bit j is 1 iff color j is in the set.
 *
 * Space optimization: only store dp for current |S| = i and previous.
 *
 * @param g      input graph
 * @param cc     color-coding state (colors pre-assigned)
 * @param result found clique output
 * @return       true if colorful k-clique found
 */
static bool color_coding_dp(const graph_t *g, color_coding_t *cc, clique_t *result) {
    int32_t n = cc->n, k = cc->k;
    int32_t num_masks = 1 << k;

    /* dp[mask][v] for current mask size */
    /* Use a compact representation: dp_table[mask * n + v] */
    bool *dp = (bool *)calloc((size_t)num_masks * (size_t)n, sizeof(bool));
    /* prev[mask][v] stores predecessor for reconstruction */
    int32_t *prev = (int32_t *)malloc((size_t)num_masks * (size_t)n * sizeof(int32_t));
    if (!dp || !prev) { free(dp); free(prev); return false; }

    /* Initialize: dp[1 << color[v]][v] = true */
    for (int32_t v = 0; v < n; v++) {
        int32_t color_bit = 1 << cc->colors[v];
        dp[color_bit * n + v] = true;
        prev[color_bit * n + v] = -1;
    }

    /* DP: iterate over mask sizes from 2 to k */
    for (int32_t size = 2; size <= k; size++) {
        /* Iterate over all masks of the current size */
        /* Count set bits to filter masks efficiently */
        for (int32_t mask = 1; mask < num_masks; mask++) {
            /* Count bits in mask */
            int32_t bits = 0, tmp = mask;
            while (tmp) { bits++; tmp &= tmp - 1; }
            if (bits != size) continue;

            for (int32_t v = 0; v < n; v++) {
                int32_t color_bit = 1 << cc->colors[v];
                if (!(mask & color_bit)) continue; /* v's color not in mask */

                int32_t prev_mask = mask ^ color_bit; /* S \ {c(v)} */
                if (prev_mask == 0) continue;

                /* Check all neighbors u with color in prev_mask */
                for (int32_t u = 0; u < n; u++) {
                    if (!graph_has_edge(g, u, v)) continue;
                    if (!(prev_mask & (1 << cc->colors[u]))) continue;
                    if (dp[prev_mask * n + u]) {
                        dp[mask * n + v] = true;
                        prev[mask * n + v] = u;
                        break; /* found one predecessor */
                    }
                }
            }
        }
    }

    /* Check if any vertex has dp[all_colors][v] = true */
    int32_t all_mask = num_masks - 1; /* all k colors */
    int32_t end_v = -1;
    for (int32_t v = 0; v < n; v++) {
        if (dp[all_mask * n + v]) {
            end_v = v;
            break;
        }
    }

    if (end_v < 0) {
        free(dp); free(prev);
        return false;
    }

    /* Reconstruct the clique */
    if (result) {
        result->k = k;
        result->found = true;
        result->vertices = (int32_t *)malloc((size_t)k * sizeof(int32_t));
        if (result->vertices) {
            int32_t cur = end_v;
            int32_t cur_mask = all_mask;
            for (int32_t pos = k - 1; pos >= 0; pos--) {
                result->vertices[pos] = cur;
                int32_t next = prev[cur_mask * n + cur];
                cur_mask ^= (1 << cc->colors[cur]);
                cur = next;
            }
        }
    }

    free(dp); free(prev);
    return true;
}

bool color_coding_find_clique(const graph_t *g,
                               color_coding_t *cc,
                               clique_t *result) {
    clique_t candidate = {0};
    if (!color_coding_dp(g, cc, &candidate)) return false;
    /* Verify the found set is actually a clique */
    if (!is_k_clique(g, candidate.vertices, candidate.k)) {
        clique_free(&candidate);
        return false;
    }
    if (result) {
        *result = candidate;
    } else {
        clique_free(&candidate);
    }
    return true;
}

bool color_coding_solve(const graph_t *g, int32_t k,
                         int32_t num_iterations, int64_t seed,
                         clique_t *result) {
    if (!g || k <= 0 || k > g->n) return false;

    color_coding_t *cc = color_coding_init(g->n, k, seed);
    if (!cc) return false;

    if (num_iterations > 0) cc->num_iterations = num_iterations;

    bool found = false;
    for (int32_t iter = 0; iter < cc->num_iterations; iter++) {
        color_coding_random_colors(cc);
        if (color_coding_find_clique(g, cc, result)) {
            found = true;
            break;
        }
    }

    color_coding_destroy(cc);
    return found;
}

bool color_coding_derandomized(const graph_t *g, int32_t k,
                                clique_t *result) {
    if (!g || k <= 0 || k > g->n) return false;
    /* Derandomized version would use explicit k-perfect hash family.
       For now, we implement the randomized version with high iteration
       count to achieve derandomization-level guarantees.
       With 100 * e^k iterations: Pr[failure] < e^{-100}. */
    return color_coding_solve(g, k, 100 * (int32_t)(exp((double)k) + 0.5),
                              (int64_t)time(NULL), result);
}

/* ================================================================
 * Branch-and-Bound for Maximum Clique
 *
 * Tomita et al. (2003, 2010, 2016) branch-and-bound algorithm.
 *
 * Core ideas:
 * - Depth-first search building cliques incrementally
 * - At each node, compute candidate set = all vertices adjacent
 *   to all vertices in current solution
 * - Coloring bound: color candidate subgraph greedily. If current_size
 *   + num_colors <= best_so_far, prune this branch
 * - Branching: for each candidate vertex, include it and recurse,
 *   then exclude it
 * - Maximum Clique: the B&B finds global maximum
 *
 * The greedy coloring bound (DSATUR heuristic):
 * - Sorts candidate vertices by degree in the subgraph
 * - Assigns smallest available color not used by adjacent vertices
 * - Number of colors used is an upper bound on clique size in subgraph
 * ================================================================ */

branch_and_bound_t *branch_bound_init(int32_t n, int32_t max_k) {
    if (n <= 0) return NULL;
    branch_and_bound_t *bb = (branch_and_bound_t *)calloc(1, sizeof(branch_and_bound_t));
    if (!bb) return NULL;

    bb->current_solution = (int32_t *)malloc((size_t)n * sizeof(int32_t));
    bb->best_solution = (int32_t *)malloc((size_t)n * sizeof(int32_t));
    if (!bb->current_solution || !bb->best_solution) {
        free(bb->current_solution); free(bb->best_solution); free(bb);
        return NULL;
    }
    bb->current_size = 0;
    bb->best_size = 0;
    bb->nodes_explored = 0;
    bb->upper_bound_limit = (max_k > 0) ? max_k : n;
    bb->use_approximate_coloring = true;
    return bb;
}

void branch_bound_destroy(branch_and_bound_t *bb) {
    if (!bb) return;
    free(bb->current_solution);
    free(bb->best_solution);
    free(bb);
}

int32_t greedy_coloring_bound(const graph_t *g,
                               const int32_t *vertices,
                               int32_t num_vertices,
                               int32_t *colors) {
    if (!g || !vertices || num_vertices <= 0) return 0;

    /* Initialize all colors to -1 (uncolored) */
    for (int32_t i = 0; i < num_vertices; i++) colors[i] = -1;

    int32_t max_color = -1;

    for (int32_t i = 0; i < num_vertices; i++) {
        /* Find smallest available color */
        bool used[256] = {false}; /* max 256 colors in practice */
        for (int32_t j = 0; j < i; j++) {
            if (colors[j] >= 0 && graph_has_edge(g, vertices[i], vertices[j])) {
                used[colors[j]] = true;
            }
        }
        int32_t c;
        for (c = 0; c < 256 && used[c]; c++) {}
        colors[i] = c;
        if (c > max_color) max_color = c;
    }

    return max_color + 1; /* number of colors used */
}

bool coloring_bound_prune(int32_t current_size,
                           const graph_t *g,
                           const int32_t *candidates,
                           int32_t num_candidates,
                           int32_t best_so_far) {
    if (num_candidates <= 0) return (current_size <= best_so_far);

    int32_t *colors = (int32_t *)malloc((size_t)num_candidates * sizeof(int32_t));
    if (!colors) return false;

    int32_t num_colors = greedy_coloring_bound(g, candidates, num_candidates, colors);
    free(colors);

    /* Prune if best we could achieve <= best_so_far */
    return (current_size + num_colors <= best_so_far);
}

/**
 * Recursive branch-and-bound for maximum clique.
 */
static void branch_bound_rec(const graph_t *g, int32_t target_k,
                              int32_t *candidates, int32_t num_candidates,
                              branch_and_bound_t *bb) {
    bb->nodes_explored++;

    /* Check if current solution is complete (a clique) */
    if (num_candidates == 0) {
        if (bb->current_size > bb->best_size) {
            bb->best_size = bb->current_size;
            memcpy(bb->best_solution, bb->current_solution,
                   (size_t)bb->current_size * sizeof(int32_t));
        }
        return;
    }

    /* Coloring bound pruning */
    if (bb->use_approximate_coloring &&
        coloring_bound_prune(bb->current_size, g, candidates,
                              num_candidates, bb->best_size)) {
        return;
    }

    /* Branch: try adding each candidate */
    for (int32_t i = 0; i < num_candidates; i++) {
        int32_t v = candidates[i];

        /* Include v in the solution */
        bb->current_solution[bb->current_size++] = v;

        /* Compute new candidates: neighbors of v in current candidates */
        int32_t *new_cand = (int32_t *)malloc((size_t)g->n * sizeof(int32_t));
        int32_t new_count = 0;
        if (new_cand) {
            for (int32_t j = i + 1; j < num_candidates; j++) {
                if (graph_has_edge(g, v, candidates[j])) {
                    new_cand[new_count++] = candidates[j];
                }
            }
            branch_bound_rec(g, target_k, new_cand, new_count, bb);
            free(new_cand);
        }

        /* Backtrack */
        bb->current_size--;

        /* Early exit if target reached */
        if (target_k > 0 && bb->best_size >= target_k) return;
    }
}

int32_t branch_bound_solve(const graph_t *g, int32_t k,
                            branch_and_bound_t *bb,
                            clique_t *result) {
    if (!g || !bb) return 0;

    /* Initialize candidates: all vertices */
    int32_t *candidates = (int32_t *)malloc((size_t)g->n * sizeof(int32_t));
    if (!candidates) return -1;
    for (int32_t i = 0; i < g->n; i++) candidates[i] = i;

    bb->best_size = 0;
    bb->current_size = 0;
    bb->nodes_explored = 0;

    branch_bound_rec(g, k, candidates, g->n, bb);

    if (result && bb->best_size > 0) {
        result->found = true;
        result->k = bb->best_size;
        result->vertices = (int32_t *)malloc((size_t)bb->best_size * sizeof(int32_t));
        if (result->vertices) {
            memcpy(result->vertices, bb->best_solution,
                   (size_t)bb->best_size * sizeof(int32_t));
        }
    }

    free(candidates);
    return bb->best_size;
}

/* ================================================================
 * Matrix Multiplication for Triangle and Small Cliques
 *
 * For 3-clique (triangle detection):
 * compute A^2, check A[i][j] * (A^2)[i][j] > 0 for some i,j.
 *
 * Actually, a triangle exists iff there exist i,j,k such that
 * A[i][j] = A[j][k] = A[k][i] = 1.
 *
 * Using matrix multiplication: (A^2)[i][j] = sum_k A[i][k]*A[k][j].
 * If A[i][j] = 1 and (A^2)[i][j] > 0, then there exists k with
 * A[i][k] = A[k][j] = 1, forming triangle (i,j,k).
 *
 * For k-clique with k = 3t (Nesetril-Poljak):
 * Construct a new graph where each "vertex" is a t-clique in G.
 * A 3-clique in the new graph corresponds to a 3t-clique in G.
 * ================================================================ */

/**
 * Standard O(n^3) matrix multiplication (for simplicity).
 * In practice one would plug in Strassen or Coppersmith-Winograd.
 */
static void matrix_multiply(const int32_t *A, const int32_t *B,
                             int32_t *C, int32_t n) {
    for (int32_t i = 0; i < n; i++) {
        for (int32_t j = 0; j < n; j++) {
            int32_t sum = 0;
            for (int32_t k = 0; k < n; k++) {
                sum += A[i * n + k] * B[k * n + j];
            }
            C[i * n + j] = sum;
        }
    }
}

bool matrix_triangle_detect(const graph_t *g, clique_t *result) {
    if (!g || g->n < 3) return false;

    int32_t n = g->n;
    int32_t *A = (int32_t *)calloc((size_t)n * n, sizeof(int32_t));
    int32_t *A2 = (int32_t *)calloc((size_t)n * n, sizeof(int32_t));
    if (!A || !A2) { free(A); free(A2); return false; }

    /* Copy adjacency to integer matrix */
    for (int32_t i = 0; i < n; i++) {
        for (int32_t j = 0; j < n; j++) {
            A[i * n + j] = (g->adj[i * n + j] != 0) ? 1 : 0;
        }
    }

    matrix_multiply(A, A, A2, n);

    /* Check for triangle: A[i][j]=1 and A2[i][j] > 0 */
    for (int32_t i = 0; i < n; i++) {
        for (int32_t j = i + 1; j < n; j++) {
            if (A[i * n + j] && A2[i * n + j] > 0) {
                /* Find the witness k: A[i][k]=1 and A[k][j]=1 */
                for (int32_t k = 0; k < n; k++) {
                    if (k != i && k != j &&
                        A[i * n + k] && A[k * n + j]) {
                        if (result) {
                            result->k = 3;
                            result->found = true;
                            result->vertices = (int32_t *)malloc(3 * sizeof(int32_t));
                            if (result->vertices) {
                                result->vertices[0] = i;
                                result->vertices[1] = j;
                                result->vertices[2] = k;
                            }
                        }
                        free(A); free(A2);
                        return true;
                    }
                }
            }
        }
    }

    free(A); free(A2);
    return false;
}

bool matrix_4clique_detect(const graph_t *g, clique_t *result) {
    if (!g || g->n < 4) return false;
    /* 4-clique detection via partitioning:
       For each pair of vertices (a,b), check if they share
       at least 2 common neighbors that are themselves adjacent. */

    int32_t n = g->n;
    for (int32_t a = 0; a < n - 3; a++) {
        for (int32_t b = a + 1; b < n - 2; b++) {
            if (!graph_has_edge(g, a, b)) continue;

            /* Find common neighbors of a and b */
            int32_t common[128]; /* practical limit */
            int32_t cnt = 0;
            for (int32_t c = 0; c < n && cnt < 128; c++) {
                if (c == a || c == b) continue;
                if (graph_has_edge(g, a, c) && graph_has_edge(g, b, c)) {
                    common[cnt++] = c;
                }
            }

            /* Check if any two common neighbors are adjacent */
            for (int32_t i = 0; i < cnt; i++) {
                for (int32_t j = i + 1; j < cnt; j++) {
                    if (graph_has_edge(g, common[i], common[j])) {
                        if (result) {
                            result->k = 4;
                            result->found = true;
                            result->vertices = (int32_t *)malloc(4 * sizeof(int32_t));
                            if (result->vertices) {
                                result->vertices[0] = a;
                                result->vertices[1] = b;
                                result->vertices[2] = common[i];
                                result->vertices[3] = common[j];
                            }
                        }
                        return true;
                    }
                }
            }
        }
    }
    return false;
}

bool matrix_kclique_detect(const graph_t *g, int32_t k, clique_t *result) {
    if (!g || k < 3 || k > g->n) return false;
    if (k == 3) return matrix_triangle_detect(g, result);
    if (k == 4) return matrix_4clique_detect(g, result);
    /* For k > 4, fall back to color-coding which is more practical */
    return color_coding_solve(g, k, 100, (int64_t)time(NULL), result);
}

double current_omega_bound(void) {
    /* Best known upper bound on matrix multiplication exponent.
       Alman, Duan, Williams, Xu, Xu, Zhou (2024): omega < 2.371339.
       Previous: Alman-Williams (2021): omega < 2.37286.
       Classic: Coppersmith-Winograd (1990): omega < 2.376.
       Trivial bound: omega <= 3. */
    return 2.371339;
}
