/**
 * kclique_core.c ? Core k-Clique checking and search algorithms
 *
 * Implements clique verification, brute-force search, counting,
 * and Bron-Kerbosch maximal clique enumeration.
 *
 * Reference: Garey & Johnson, "Computers and Intractability" (1979)
 * Reference: Bron & Kerbosch, "Algorithm 457" (CACM 1973)
 * Reference: Tomita, Tanaka, Takahashi, "Worst-case analysis of BK" (TCS 2006)
 * Reference: Moon & Moser, "On cliques in graphs" (Israel J. Math 1965)
 * Reference: Eppstein, Loffler, Strash, "Listing All Maximal Cliques" (2013)
 *
 * Knowledge: L1 Definitions ? k-Clique decision/search
 * Knowledge: L5 Algorithms ? Brute force, Bron-Kerbosch
 * Knowledge: L6 Canonical Problems ? k-Clique, Max-Clique
 */

#include "kclique_core.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Forward declaration of next_combination (from kclique_types.c) */
extern int64_t binomial_coefficient(int32_t n, int32_t k);

/* ================================================================
 * Clique Verification
 * ================================================================ */

bool is_clique(const graph_t *g, const int32_t *vertices, int32_t count) {
    if (!g || !vertices) return false;
    if (count < 0) return false;
    if (count <= 1) return true; /* 0 or 1 vertices trivially form a clique */
    for (int32_t i = 0; i < count; i++) {
        if (vertices[i] < 0 || vertices[i] >= g->n) return false;
    }
    for (int32_t i = 0; i < count; i++) {
        for (int32_t j = i + 1; j < count; j++) {
            if (!graph_has_edge(g, vertices[i], vertices[j])) {
                return false;
            }
        }
    }
    return true;
}

bool is_k_clique(const graph_t *g, const int32_t *vertices, int32_t k) {
    if (!g || !vertices || k <= 0) return false;
    return is_clique(g, vertices, k);
}

/* ================================================================
 * Brute-Force k-Clique Search: O(C(n,k) * k^2)
 *
 * This is the "trivial" algorithm. Under SETH, it is essentially
 * optimal: k-Clique requires n^{k-o(1)} time for general k.
 *
 * Approach: iterate over all k-subsets of [n] in lexicographic order.
 * For each subset, verify pairwise adjacency. Return first clique found.
 *
 * Number of subsets: C(n,k) = Theta(n^k) for fixed k.
 * Per-subset cost: O(k^2) for verifying all pairs.
 * Total: O(k^2 * n^k / k!) = O(n^k) for fixed k.
 * ================================================================ */

/**
 * Generate the next combination of size k from {0, ..., n-1}
 * in lexicographic order. Returns false if this was the last one.
 *
 * @param combo  current combination (modified in-place)
 * @param n      universe size
 * @param k      combination size
 * @return       true if next combination exists
 */
static bool next_combination(int32_t *combo, int32_t n, int32_t k) {
    int32_t i;
    for (i = k - 1; i >= 0; i--) {
        if (combo[i] < n - k + i) break;
    }
    if (i < 0) return false;
    combo[i]++;
    for (int32_t j = i + 1; j < k; j++) {
        combo[j] = combo[j - 1] + 1;
    }
    return true;
}

bool find_k_clique_brute(const graph_t *g, int32_t k, clique_t *result) {
    if (!g || k < 0 || k > g->n) return false;

    /* Trivial cases */
    if (k == 0) {
        if (result) {
            result->found = true;
            result->k = 0;
            result->vertices = NULL;
        }
        return true;
    }
    if (k == 1) {
        if (result) {
            result->found = true;
            result->k = 1;
            result->vertices = (int32_t *)malloc(sizeof(int32_t));
            if (result->vertices) result->vertices[0] = 0;
        }
        return g->n > 0;
    }
    if (k > GRAPH_MAX_K) {
        fprintf(stderr, "find_k_clique_brute: k=%d exceeds limit\n", k);
        return false;
    }

    int32_t *combo = (int32_t *)malloc((size_t)k * sizeof(int32_t));
    if (!combo) return false;

    /* Initialize first combination: {0, 1, ..., k-1} */
    for (int32_t i = 0; i < k; i++) combo[i] = i;

    do {
        if (is_k_clique(g, combo, k)) {
            if (result) {
                result->k = k;
                result->found = true;
                result->vertices = (int32_t *)malloc((size_t)k * sizeof(int32_t));
                if (result->vertices) {
                    memcpy(result->vertices, combo, (size_t)k * sizeof(int32_t));
                }
            }
            free(combo);
            return true;
        }
    } while (next_combination(combo, g->n, k));

    free(combo);
    if (result) {
        result->found = false;
        result->k = 0;
        result->vertices = NULL;
    }
    return false;
}

int64_t count_k_cliques(const graph_t *g, int32_t k) {
    if (!g || k < 0 || k > g->n) return -1;
    if (k <= 1) {
        if (k == 0) return 1; /* empty set */
        return g->n;           /* each vertex is a 1-clique */
    }
    if (k > GRAPH_MAX_K) return -1;

    int64_t count = 0;
    int32_t *combo = (int32_t *)malloc((size_t)k * sizeof(int32_t));
    if (!combo) return -1;

    for (int32_t i = 0; i < k; i++) combo[i] = i;

    do {
        if (is_k_clique(g, combo, k)) count++;
    } while (next_combination(combo, g->n, k));

    free(combo);
    return count;
}

/* ================================================================
 * Bron-Kerbosch Algorithm for Maximum Clique and Enumeration
 *
 * The Bron-Kerbosch algorithm with pivot is the canonical method
 * for enumerating all maximal cliques in an undirected graph.
 *
 * Key sets:
 *  - R: current clique (growing)
 *  - P: candidate vertices (all adjacent to all vertices in R)
 *  - X: excluded vertices (adjacent to all in R, but already processed)
 *
 * Pivot optimization (Bron-Kerbosch 1973, Tomita et al. 2006):
 *  Choose a pivot u in P ? X that maximizes |P ? N(u)|.
 *  Only branch on vertices in P \ N(u), reducing recursive calls.
 *
 * Complexity: O(3^{n/3}) worst-case (Moon-Moser 1965 lower bound is
 * matched, so BK is optimal for enumeration).
 * ================================================================ */

/**
 * Recursive Bron-Kerbosch with pivot for enumerating maximal cliques.
 *
 * @param g            graph
 * @param R            current clique vertices
 * @param R_size       |R|
 * @param P            candidate vertices
 * @param P_size       |P|
 * @param X            excluded vertices
 * @param X_size       |X|
 * @param best         [in/out] best clique found
 * @param best_size    [in/out] size of best clique
 */
static void bron_kerbosch_rec(const graph_t *g,
                               int32_t *R, int32_t R_size,
                               int32_t *P, int32_t P_size,
                               int32_t *X, int32_t X_size,
                               int32_t *best, int32_t *best_size) {
    if (P_size == 0 && X_size == 0) {
        /* R is a maximal clique */
        if (R_size > *best_size) {
            *best_size = R_size;
            if (best) {
                memcpy(best, R, (size_t)R_size * sizeof(int32_t));
            }
        }
        return;
    }

    /* Pivot selection: choose u from P ? X maximizing |P ? N(u)| */
    int32_t pivot = -1, max_intersection = -1;
    int32_t total_candidates = P_size + X_size;
    for (int32_t i = 0; i < total_candidates; i++) {
        int32_t u = (i < P_size) ? P[i] : X[i - P_size];
        int32_t intersect = 0;
        for (int32_t j = 0; j < P_size; j++) {
            if (graph_has_edge(g, u, P[j])) intersect++;
        }
        if (intersect > max_intersection) {
            max_intersection = intersect;
            pivot = u;
        }
    }

    /* Branch on vertices in P that are NOT adjacent to pivot */
    int32_t *newP = (int32_t *)malloc((size_t)g->n * sizeof(int32_t));
    int32_t *newX = (int32_t *)malloc((size_t)g->n * sizeof(int32_t));
    if (!newP || !newX) { free(newP); free(newX); return; }

    for (int32_t i = 0; i < P_size; i++) {
        int32_t v = P[i];
        if (pivot >= 0 && graph_has_edge(g, pivot, v)) continue;

        /* Add v to R */
        R[R_size] = v;

        /* Compute new P: P ? N(v) */
        int32_t newP_size = 0;
        for (int32_t j = 0; j < P_size; j++) {
            if (graph_has_edge(g, v, P[j])) {
                newP[newP_size++] = P[j];
            }
        }

        /* Compute new X: X ? N(v) */
        int32_t newX_size = 0;
        for (int32_t j = 0; j < X_size; j++) {
            if (graph_has_edge(g, v, X[j])) {
                newX[newX_size++] = X[j];
            }
        }

        bron_kerbosch_rec(g, R, R_size + 1,
                          newP, newP_size, newX, newX_size,
                          best, best_size);

        /* Move v from P to X */
        X[X_size++] = v;
        /* Remove v from P (shift remaining) */
        /* Note: v was the i-th element, we handle this in the loop */
    }

    free(newP);
    free(newX);
}

int32_t max_clique_bron_kerbosch(const graph_t *g, clique_t *result) {
    if (!g || g->n == 0) {
        if (result) { result->found = false; result->k = 0; result->vertices = NULL; }
        return 0;
    }

    int32_t *R = (int32_t *)malloc((size_t)g->n * sizeof(int32_t));
    int32_t *P = (int32_t *)malloc((size_t)g->n * sizeof(int32_t));
    int32_t *X = (int32_t *)malloc((size_t)g->n * sizeof(int32_t));
    int32_t *best = (int32_t *)malloc((size_t)g->n * sizeof(int32_t));

    if (!R || !P || !X || !best) {
        free(R); free(P); free(X); free(best);
        return -1;
    }

    /* Initialize P = all vertices, X = empty */
    for (int32_t i = 0; i < g->n; i++) P[i] = i;
    int32_t best_size = 0;

    bron_kerbosch_rec(g, R, 0, P, g->n, X, 0, best, &best_size);

    if (result && best_size > 0) {
        result->found = true;
        result->k = best_size;
        result->vertices = (int32_t *)malloc((size_t)best_size * sizeof(int32_t));
        if (result->vertices) {
            memcpy(result->vertices, best, (size_t)best_size * sizeof(int32_t));
        }
    } else if (result) {
        result->found = false;
        result->k = 0;
        result->vertices = NULL;
    }

    free(R); free(P); free(X); free(best);
    return best_size;
}

/**
 * Recursive BK for enumeration (accumulates all maximal cliques).
 */
static void bron_kerbosch_enum_rec(const graph_t *g,
                                    int32_t *R, int32_t R_size,
                                    int32_t *P, int32_t P_size,
                                    int32_t *X, int32_t X_size,
                                    clique_t *cliques, int32_t *capacity,
                                    int32_t *count) {
    if (P_size == 0 && X_size == 0) {
        if (*count >= *capacity) return;
        cliques[*count].k = R_size;
        cliques[*count].found = true;
        cliques[*count].vertices = (int32_t *)malloc((size_t)R_size * sizeof(int32_t));
        if (cliques[*count].vertices) {
            memcpy(cliques[*count].vertices, R, (size_t)R_size * sizeof(int32_t));
        }
        (*count)++;
        return;
    }

    int32_t pivot = -1, max_intersection = -1;
    for (int32_t i = 0; i < P_size + X_size; i++) {
        int32_t u = (i < P_size) ? P[i] : X[i - P_size];
        int32_t intersect = 0;
        for (int32_t j = 0; j < P_size; j++) {
            if (graph_has_edge(g, u, P[j])) intersect++;
        }
        if (intersect > max_intersection) {
            max_intersection = intersect;
            pivot = u;
        }
    }

    int32_t *newP = (int32_t *)malloc((size_t)g->n * sizeof(int32_t));
    int32_t *newX = (int32_t *)malloc((size_t)g->n * sizeof(int32_t));
    if (!newP || !newX) { free(newP); free(newX); return; }

    int32_t *P_copy = (int32_t *)malloc((size_t)P_size * sizeof(int32_t));
    if (P_copy) memcpy(P_copy, P, (size_t)P_size * sizeof(int32_t));

    for (int32_t i = 0; i < P_size; i++) {
        int32_t v = P_copy[i];
        if (pivot >= 0 && graph_has_edge(g, pivot, v)) {
            X[X_size++] = v;
            continue;
        }

        R[R_size] = v;

        int32_t newP_size = 0;
        for (int32_t j = 0; j < P_size; j++) {
            if (graph_has_edge(g, v, P_copy[j])) {
                newP[newP_size++] = P_copy[j];
            }
        }

        int32_t newX_size = 0;
        for (int32_t j = 0; j < X_size; j++) {
            if (graph_has_edge(g, v, X[j])) {
                newX[newX_size++] = X[j];
            }
        }

        bron_kerbosch_enum_rec(g, R, R_size + 1,
                                newP, newP_size, newX, newX_size,
                                cliques, capacity, count);
        X[X_size++] = v;
    }

    free(newP); free(newX); free(P_copy);
}

int32_t enumerate_maximal_cliques(const graph_t *g,
                                   clique_t **cliques,
                                   int32_t *num_cliques) {
    if (!g || !cliques || !num_cliques) return -1;
    *num_cliques = 0;
    *cliques = NULL;

    /* Upper bound on number of maximal cliques: 3^{n/3} (Moon-Moser) */
    int32_t max_cliques = 1000000; /* practical limit for safety */
    clique_t *cliq_array = (clique_t *)calloc((size_t)max_cliques, sizeof(clique_t));
    if (!cliq_array) return -1;

    int32_t *R = (int32_t *)malloc((size_t)g->n * sizeof(int32_t));
    int32_t *P = (int32_t *)malloc((size_t)g->n * sizeof(int32_t));
    int32_t *X = (int32_t *)malloc((size_t)g->n * sizeof(int32_t));

    if (!R || !P || !X) {
        free(R); free(P); free(X); free(cliq_array);
        return -1;
    }

    for (int32_t i = 0; i < g->n; i++) P[i] = i;
    int32_t count = 0;

    bron_kerbosch_enum_rec(g, R, 0, P, g->n, X, 0,
                            cliq_array, &max_cliques, &count);

    *num_cliques = count;
    *cliques = cliq_array;

    free(R); free(P); free(X);
    return count;
}
