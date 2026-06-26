/**
 * kclique_count.c ? Advanced clique counting algorithms
 *
 * Implements degeneracy-based counting, neighborhood intersection,
 * inclusion-exclusion, and other techniques for counting cliques
 * beyond the naive O(n^k) enumeration.
 *
 * Reference: Chiba & Nishizeki, "Arboricity and Subgraph Listing" (SICOMP 1985)
 * Reference: Eppstein, "Subgraph Isomorphism in Planar Graphs" (JGT 1999)
 * Reference: Alon, Dao, Hajirasouliha, Hormozdiari, Sahinalp,
 *             "Biomolecular network motif counting" (Bioinformatics 2008)
 * Reference: Jain & Seshadhri, "A Fast and Provable Method for Estimating
 *             Clique Counts" (KDD 2017)
 *
 * Knowledge: L5 Algorithms ? Advanced clique counting methods
 * Knowledge: L7 Applications ? Motif counting in biological networks
 * Knowledge: L8 Advanced Topics ? Sublinear-space counting
 */

#include "kclique_types.h"
#include "kclique_core.h"
#include "kclique_graph.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ================================================================
 * Degeneracy-Based Clique Counting
 *
 * Key insight: for each vertex v, the number of k-cliques containing
 * v equals the number of (k-1)-cliques in the induced subgraph of N(v).
 *
 * By processing vertices in a degeneracy ordering, we can bound the
 * size of N(v) by the degeneracy d (typically much smaller than n).
 *
 * Algorithm (Chiba-Nishizeki 1985):
 * 1. Compute degeneracy ordering v_1, ..., v_n
 * 2. For each v_i in order:
 *    a. Let S = N(v_i) ? {v_{i+1}, ..., v_n}
 *    b. Count (k-1)-cliques in the subgraph induced by S
 *    c. Each such (k-1)-clique together with v_i forms a k-clique
 *
 * Complexity: O(k * d^{k-2} * m) where d = degeneracy.
 * This beats O(n^k) when d << n (which is common in sparse graphs).
 * ================================================================ */

/**
 * Count (k-1)-cliques in a subgraph induced by a vertex set using
 * a simple recursive approach for small k.
 *
 * @param g          original graph
 * @param vertices   vertex subset to search in
 * @param count      size of the vertex subset
 * @param k_minus_1  target clique size (k-1)
 * @return           number of (k-1)-cliques in the induced subgraph
 */
static int64_t count_cliques_in_subset(const graph_t *g,
                                        const int32_t *vertices,
                                        int32_t count,
                                        int32_t k_minus_1) {
    if (k_minus_1 <= 0) return 1; /* empty clique counted once */
    if (k_minus_1 == 1) return count; /* each vertex is a 1-clique */
    if (count < k_minus_1) return 0;

    int64_t total = 0;
    /* Recursively count: for each vertex, count cliques containing it
       in its neighborhood within the remaining set. */
    for (int32_t i = 0; i <= count - k_minus_1; i++) {
        int32_t v = vertices[i];
        /* Find neighbors of v among later vertices */
        int32_t later_neighbors[256]; /* practical limit for small k */
        int32_t ln_count = 0;
        for (int32_t j = i + 1; j < count && ln_count < 256; j++) {
            if (graph_has_edge(g, v, vertices[j])) {
                later_neighbors[ln_count++] = vertices[j];
            }
        }
        total += count_cliques_in_subset(g, later_neighbors, ln_count,
                                          k_minus_1 - 1);
    }
    return total;
}

/**
 * Count k-cliques using the degeneracy ordering method.
 *
 * @param g  graph
 * @param k  clique size
 * @return   number of k-cliques, or -1 on error
 */
int64_t count_k_cliques_degeneracy(const graph_t *g, int32_t k) {
    if (!g || k < 0 || k > g->n) return -1;
    if (k <= 1) return (k == 0) ? 1 : g->n;

    /* Compute degeneracy ordering */
    int32_t *ordering;
    int32_t degeneracy_val = graph_degeneracy(g, &ordering);
    if (degeneracy_val < 0 || !ordering) return -1;

    int64_t total = 0;
    int32_t *later_set = (int32_t *)malloc((size_t)g->n * sizeof(int32_t));
    if (!later_set) { free(ordering); return -1; }

    for (int32_t pos = 0; pos < g->n; pos++) {
        int32_t v = ordering[pos];
        /* Collect vertices later in the ordering that are adjacent to v */
        int32_t later_count = 0;
        for (int32_t p = pos + 1; p < g->n; p++) {
            if (graph_has_edge(g, v, ordering[p])) {
                later_set[later_count++] = ordering[p];
            }
        }
        /* Count (k-1)-cliques in the later neighborhood */
        total += count_cliques_in_subset(g, later_set, later_count, k - 1);
    }

    free(later_set);
    free(ordering);
    return total;
}

/* ================================================================
 * Neighborhood Intersection Method
 *
 * For counting k-cliques, we can use:
 *
 * count = sum_{v_1} ... sum_{v_k} product_{i<j} I[(v_i, v_j) ? E] / k!
 *
 * Optimized: for each edge (u,v), count triangles sharing this edge,
 * then for each triangle count 4-cliques containing it, etc.
 *
 * This is the "edge-iterator" method: iterate over cliques of size s,
 * then extend to size s+1 by checking common neighbors.
 * ================================================================ */

/**
 * Count k-cliques using the Chiba-Nishizeki edge-iteration method.
 *
 * For k=3: iterate edges, for each edge (u,v), count common neighbors.
 * For k=4: iterate triangles, for each triangle, count common neighbors.
 * For general k: iterate (k-1)-cliques, for each, count common neighbors.
 *
 * Complexity: O(k * m * d^{k-2}) where d = degeneracy.
 *
 * @param g  graph
 * @param k  clique size
 * @return   number of k-cliques
 */
int64_t count_k_cliques_iterative(const graph_t *g, int32_t k) {
    if (!g || k < 0 || k > g->n) return -1;
    if (k <= 1) return (k == 0) ? 1 : g->n;
    if (k == 2) return g->m; /* each edge is a 2-clique */
    if (k == 3) {
        /* Count triangles: sum over edges of |N(u) ? N(v)| */
        int64_t triangles = 0;
        for (int32_t u = 0; u < g->n; u++) {
            for (int32_t v = u + 1; v < g->n; v++) {
                if (!graph_has_edge(g, u, v)) continue;
                /* Count common neighbors of u and v */
                for (int32_t w = v + 1; w < g->n; w++) {
                    if (graph_has_edge(g, u, w) && graph_has_edge(g, v, w)) {
                        triangles++;
                    }
                }
            }
        }
        return triangles;
    }

    /* For k >= 4, use recursive extension from smaller cliques.
       This is essentially the same as the subset counting but with
       a different iteration order (edge-first rather than vertex-first). */
    return count_k_cliques_degeneracy(g, k);
}

/* ================================================================
 * Inclusion-Exclusion Counting
 *
 * For a subset S of vertices, define f(S) = number of cliques
 * contained entirely within S. By inclusion-exclusion:
 *
 * Number of k-cliques = sum_{T subset of [n], |T| >= k}
 *                       (-1)^{|T|-k} * C(|T|-1, k-1) * f(T)
 *
 * This is generally impractical but provides a theoretical framework.
 * ================================================================ */

/**
 * Compute the number of k-cliques using explicit inclusion-exclusion.
 * This is for theoretical illustration; O(2^n) complexity makes it
 * impractical for n > 30.
 *
 * @param g  graph (n <= 20 recommended)
 * @param k  clique size
 * @return   number of k-cliques, or -1 on error
 */
int64_t count_k_cliques_inclusion_exclusion(const graph_t *g, int32_t k) {
    if (!g || k < 0 || k > g->n) return -1;
    if (g->n > 20) {
        fprintf(stderr, "Inclusion-exclusion impractical for n=%d\n", g->n);
        return -1;
    }
    if (k <= 1) return (k == 0) ? 1 : g->n;

    int64_t total = 0;
    int32_t max_mask = 1 << g->n;

    for (int32_t mask = 0; mask < max_mask; mask++) {
        /* Count bits in mask */
        int32_t bits = 0;
        int32_t tmp = mask;
        while (tmp) { bits++; tmp &= tmp - 1; }

        if (bits < k) continue;

        /* Collect vertices in this subset */
        int32_t subset[32]; /* n <= 20, so 32 is safe */
        int32_t cnt = 0;
        for (int32_t i = 0; i < g->n; i++) {
            if (mask & (1 << i)) subset[cnt++] = i;
        }

        /* Count cliques in this subset: check all k-subsets */
        int64_t cliques_in_subset = 0;
        int32_t *combo = (int32_t *)malloc((size_t)k * sizeof(int32_t));
        if (!combo) continue;
        for (int32_t i = 0; i < k; i++) combo[i] = i;

        /* Simple: check all C(cnt,k) subsets using nested loops for small k */
        if (cnt >= k) {
            /* Use brute force for subset enumeration */
            for (int32_t i1 = 0; i1 < cnt - (k-1); i1++) {
                for (int32_t i2 = i1 + 1; i2 < cnt - (k-2); i2++) {
                    if (k == 2) {
                        if (graph_has_edge(g, subset[i1], subset[i2]))
                            cliques_in_subset++;
                    } else if (k == 3) {
                        for (int32_t i3 = i2 + 1; i3 < cnt; i3++) {
                            if (graph_has_edge(g, subset[i1], subset[i2]) &&
                                graph_has_edge(g, subset[i1], subset[i3]) &&
                                graph_has_edge(g, subset[i2], subset[i3]))
                                cliques_in_subset++;
                        }
                    } else if (k == 4) {
                        for (int32_t i3 = i2 + 1; i3 < cnt - 0; i3++) {
                            for (int32_t i4 = i3 + 1; i4 < cnt; i4++) {
                                if (graph_has_edge(g, subset[i1], subset[i2]) &&
                                    graph_has_edge(g, subset[i1], subset[i3]) &&
                                    graph_has_edge(g, subset[i1], subset[i4]) &&
                                    graph_has_edge(g, subset[i2], subset[i3]) &&
                                    graph_has_edge(g, subset[i2], subset[i4]) &&
                                    graph_has_edge(g, subset[i3], subset[i4]))
                                    cliques_in_subset++;
                            }
                        }
                    }
                }
            }
        }
        free(combo);

        /* Inclusion-exclusion sign: (-1)^{|T|-k} */
        int32_t sign = ((bits - k) % 2 == 0) ? 1 : -1;
        total += sign * cliques_in_subset;
    }

    return total;
}

/**
 * Estimate the number of k-cliques using sampling.
 *
 * Randomly sample k-subsets and extrapolate from the proportion
 * that form cliques.
 *
 * @param g       graph
 * @param k       clique size
 * @param samples number of samples
 * @return        estimated count, or -1 on error
 */
double estimate_k_cliques_sampling(const graph_t *g, int32_t k, int32_t samples) {
    if (!g || k <= 0 || k > g->n || samples <= 0) return -1.0;

    int64_t hits = 0;
    int32_t *combo = (int32_t *)malloc((size_t)k * sizeof(int32_t));
    if (!combo) return -1.0;

    /* Use a simple LCG for sampling */
    uint64_t seed = 12345ULL;

    for (int32_t s = 0; s < samples; s++) {
        /* Reservoir-like sampling of k distinct vertices */
        bool used[4096] = {false}; /* works for n <= 4096 */
        for (int32_t i = 0; i < k; i++) {
            int32_t v;
            do {
                seed = seed * 1103515245ULL + 12345ULL;
                v = (int32_t)(seed % (uint64_t)g->n);
            } while (used[v]);
            used[v] = true;
            combo[i] = v;
        }

        if (is_k_clique(g, combo, k)) hits++;
    }

    free(combo);

    /* Extrapolate: hits/samples * C(n,k) */
    int64_t n_choose_k = 1;
    for (int32_t i = 1; i <= k; i++) {
        n_choose_k = n_choose_k * (g->n - k + i) / i;
    }

    return (double)hits * (double)n_choose_k / (double)samples;
}
