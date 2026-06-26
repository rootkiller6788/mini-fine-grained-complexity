/**
 * kclique_graph.c ? Graph generation and property analysis
 *
 * Implements Erdos-Renyi random graphs, structured graph generation,
 * degeneracy computation, core decomposition, and Turan bounds.
 *
 * Reference: Erdos & Renyi, "On Random Graphs I" (1959)
 * Reference: Bollobas, "Random Graphs" (2001)
 * Reference: Turan, "On an extremal problem in graph theory" (1941)
 * Reference: Alon & Spencer, "The Probabilistic Method" (2016)
 *
 * Knowledge: L3 Mathematical Structures ? Graph classes and properties
 * Knowledge: L5 Algorithms ? Graph generation, degeneracy ordering
 * Knowledge: L7 Applications ? Random graph models for benchmarking
 */

#include "kclique_graph.h"
#include "kclique_core.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

/* LCG for reproducibility (mirrors kclique_algorithm.c) */
static uint64_t rng_state = 123456789ULL;

static void rng_seed(uint64_t seed) {
    rng_state = seed ? seed : (uint64_t)time(NULL);
}

static uint64_t rng_next(void) {
    rng_state = rng_state * 6364136223846793005ULL + 1442695040888963407ULL;
    return rng_state;
}

static double rng_double(void) {
    return (double)(rng_next() >> 11) / (double)(1ULL << 53);
}

static int32_t __attribute__((unused)) rng_int(int32_t max) {
    if (max <= 0) return 0;
    return (int32_t)(rng_next() % (uint64_t)max);
}

/* ================================================================
 * Random Graph Generation: Erdos-Renyi G(n,p)
 *
 * In G(n,p), each of the n*(n-1)/2 possible undirected edges
 * is present independently with probability p.
 *
 * Expected edges: p * n*(n-1)/2
 * Expected degree: p*(n-1)
 * Clique number (p constant): ~ 2*log_{1/p}(n) (Bollobas 2001)
 * For p = 1/2: omega(G) ~ 2*log_2(n) asymptotically almost surely.
 * ================================================================ */

graph_t *graph_generate_erdos_renyi(const er_graph_params_t *params) {
    if (!params || params->n <= 0 || params->n > GRAPH_MAX_VERTICES) return NULL;
    if (params->p < 0.0 || params->p > 1.0) return NULL;

    rng_seed((uint64_t)params->seed);
    graph_t *g = graph_create(params->n);
    if (!g) return NULL;

    for (int32_t i = 0; i < params->n; i++) {
        for (int32_t j = i + 1; j < params->n; j++) {
            if (rng_double() < params->p) {
                graph_add_edge(g, i, j);
            }
        }
    }

    /* Plant a hidden clique if requested */
    if (params->include_clique && params->planted_k > 0 &&
        params->planted_k <= params->n) {
        /* Plant a k-clique on the first planted_k vertices */
        for (int32_t i = 0; i < params->planted_k; i++) {
            for (int32_t j = i + 1; j < params->planted_k; j++) {
                graph_add_edge(g, i, j);
            }
        }
    }

    /* Force connectivity by adding a spanning tree if needed */
    if (params->force_connected) {
        for (int32_t i = 1; i < params->n; i++) {
            if (!graph_has_edge(g, i - 1, i)) {
                graph_add_edge(g, i - 1, i);
            }
        }
    }

    return g;
}

graph_t *graph_generate_complete(int32_t n) {
    if (n <= 0 || n > GRAPH_MAX_VERTICES) return NULL;
    graph_t *g = graph_create(n);
    if (!g) return NULL;
    for (int32_t i = 0; i < n; i++) {
        for (int32_t j = i + 1; j < n; j++) {
            graph_add_edge(g, i, j);
        }
    }
    return g;
}

graph_t *graph_generate_empty(int32_t n) {
    if (n <= 0 || n > GRAPH_MAX_VERTICES) return NULL;
    return graph_create(n);
}

graph_t *graph_generate_planted_clique(int32_t n, int32_t k,
                                        double p_intra, double p_inter,
                                        int64_t seed) {
    if (n <= 0 || k <= 0 || k > n) return NULL;

    rng_seed((uint64_t)seed);
    graph_t *g = graph_create(n);
    if (!g) return NULL;

    /* Partition vertices into k groups of size n/k */
    int32_t group_size = n / k;

    for (int32_t i = 0; i < n; i++) {
        int32_t gi = i / group_size;
        if (gi >= k) gi = k - 1; /* last group may be larger */
        for (int32_t j = i + 1; j < n; j++) {
            int32_t gj = j / group_size;
            if (gj >= k) gj = k - 1;
            double p = (gi == gj) ? p_intra : p_inter;
            if (rng_double() < p) {
                graph_add_edge(g, i, j);
            }
        }
    }

    /* Plant a k-clique: one vertex from each group */
    for (int32_t gi = 0; gi < k; gi++) {
        int32_t vi = gi * group_size; /* first vertex in each group */
        for (int32_t gj = gi + 1; gj < k; gj++) {
            int32_t vj = gj * group_size;
            if (!graph_has_edge(g, vi, vj)) {
                graph_add_edge(g, vi, vj);
            }
        }
    }

    return g;
}

graph_t *graph_generate_circulant(int32_t n,
                                   const int32_t *jumps,
                                   int32_t num_jumps) {
    if (n <= 0 || n > GRAPH_MAX_VERTICES || !jumps || num_jumps <= 0) return NULL;
    graph_t *g = graph_create(n);
    if (!g) return NULL;

    for (int32_t i = 0; i < n; i++) {
        for (int32_t j_idx = 0; j_idx < num_jumps; j_idx++) {
            int32_t j = (i + jumps[j_idx]) % n;
            if (j > i) { /* add each edge once */
                graph_add_edge(g, i, j);
            }
        }
    }
    return g;
}

graph_t *graph_generate_complement(const graph_t *g) {
    return graph_complement(g);
}

/* ================================================================
 * Graph Properties
 * ================================================================ */

int32_t graph_degeneracy(const graph_t *g, int32_t **ordering) {
    if (!g) return -1;

    int32_t n = g->n;
    int32_t *deg = (int32_t *)malloc((size_t)n * sizeof(int32_t));
    int32_t *order = (int32_t *)malloc((size_t)n * sizeof(int32_t));
    bool *removed = (bool *)calloc((size_t)n, sizeof(bool));
    if (!deg || !order || !removed) {
        free(deg); free(order); free(removed);
        return -1;
    }

    /* Compute initial degrees */
    for (int32_t v = 0; v < n; v++) {
        deg[v] = graph_degree(g, v);
        if (deg[v] < 0) deg[v] = 0;
    }

    int32_t max_deg = 0;
    for (int32_t pos = 0; pos < n; pos++) {
        /* Find vertex with minimum degree among remaining */
        int32_t min_deg = n + 1, min_v = -1;
        for (int32_t v = 0; v < n; v++) {
            if (!removed[v] && deg[v] < min_deg) {
                min_deg = deg[v];
                min_v = v;
            }
        }
        if (min_v < 0) break;

        if (min_deg > max_deg) max_deg = min_deg;
        order[pos] = min_v;
        removed[min_v] = true;

        /* Update degrees of neighbors */
        for (int32_t u = 0; u < n; u++) {
            if (!removed[u] && graph_has_edge(g, min_v, u)) {
                deg[u]--;
            }
        }
    }

    if (ordering) *ordering = order;
    else free(order);
    free(deg); free(removed);
    return max_deg;
}

void graph_core_decomposition(const graph_t *g, int32_t *core_values) {
    if (!g || !core_values) return;

    int32_t n = g->n;
    int32_t *deg = (int32_t *)malloc((size_t)n * sizeof(int32_t));
    bool *removed = (bool *)calloc((size_t)n, sizeof(bool));
    if (!deg || !removed) { free(deg); free(removed); return; }

    for (int32_t v = 0; v < n; v++) {
        core_values[v] = 0;
        deg[v] = graph_degree(g, v);
    }

    /* Matula & Beck (1983) algorithm: repeatedly remove min-degree vertex */
    for (int32_t pass = 0; pass < n; pass++) {
        int32_t min_deg = n + 1, min_v = -1;
        for (int32_t v = 0; v < n; v++) {
            if (!removed[v] && deg[v] < min_deg) {
                min_deg = deg[v];
                min_v = v;
            }
        }
        if (min_v < 0) break;

        core_values[min_v] = min_deg;
        removed[min_v] = true;

        for (int32_t u = 0; u < n; u++) {
            if (!removed[u] && graph_has_edge(g, min_v, u)) {
                deg[u]--;
            }
        }
    }

    free(deg); free(removed);
}

void graph_degree_distribution(const graph_t *g,
                                int32_t *deg_count,
                                int32_t *max_deg) {
    if (!g || !deg_count) return;

    int32_t max_d = 0;
    for (int32_t v = 0; v < g->n; v++) {
        int32_t d = graph_degree(g, v);
        if (d < 0) d = 0;
        if (d < g->n) deg_count[d]++;
        if (d > max_d) max_d = d;
    }
    if (max_deg) *max_deg = max_d;
}

double graph_clustering_coefficient(const graph_t *g, int32_t v) {
    if (!g || v < 0 || v >= g->n) return -1.0;

    int32_t ncnt;
    int32_t *neighbors = graph_neighbors(g, v, &ncnt);
    if (!neighbors || ncnt < 2) {
        free(neighbors);
        return 0.0;
    }

    int32_t edges_between_neighbors = 0;
    for (int32_t i = 0; i < ncnt; i++) {
        for (int32_t j = i + 1; j < ncnt; j++) {
            if (graph_has_edge(g, neighbors[i], neighbors[j])) {
                edges_between_neighbors++;
            }
        }
    }

    free(neighbors);
    double denom = (double)ncnt * (ncnt - 1) / 2.0;
    return (denom > 0) ? (double)edges_between_neighbors / denom : 0.0;
}

double graph_average_clustering(const graph_t *g) {
    if (!g || g->n == 0) return 0.0;
    double sum = 0.0;
    for (int32_t v = 0; v < g->n; v++) {
        sum += graph_clustering_coefficient(g, v);
    }
    return sum / (double)g->n;
}

bool graph_is_triangle_free(const graph_t *g) {
    if (!g) return true;
    if (g->n < 3) return true;
    /* Check all triples */
    for (int32_t i = 0; i < g->n - 2; i++) {
        for (int32_t j = i + 1; j < g->n - 1; j++) {
            if (!graph_has_edge(g, i, j)) continue;
            for (int32_t k = j + 1; k < g->n; k++) {
                if (graph_has_edge(g, i, k) && graph_has_edge(g, j, k)) {
                    return false;
                }
            }
        }
    }
    return true;
}

int64_t turan_bound_edges(int32_t n, int32_t r) {
    if (n <= 0 || r <= 0) return 0;
    /* Turan's theorem: ex(n, K_{r+1}) = (1 - 1/r) * n^2/2 */
    /* The extremal graph is the Turan graph T(n,r): complete r-partite
       with parts as equal as possible. */
    double edges = (1.0 - 1.0 / (double)r) * (double)n * (double)(n - 1) / 2.0;
    return (int64_t)(edges + 0.5);
}

bool turan_guarantees_clique(const graph_t *g, int32_t k) {
    if (!g || k <= 1) return true;
    /* Turan's theorem: if |E| > turan_bound(n, k-1),
       then g MUST contain K_k. */
    int64_t t_bound = turan_bound_edges(g->n, k - 1);
    return (int64_t)g->m > t_bound;
}

double graph_density(const graph_t *g) {
    if (!g || g->n <= 1) return 0.0;
    double max_edges = (double)g->n * (double)(g->n - 1) / 2.0;
    return (double)g->m / max_edges;
}
