/**
 * kclique_types.c ? Graph operations and configuration initialization
 *
 * Implementation of graph creation, destruction, cloning,
 * and configuration management for the k-Clique module.
 *
 * Reference: Downey & Fellows, "Parameterized Complexity" (1999)
 * Reference: Sipser, "Introduction to the Theory of Computation" (2013)
 *
 * Knowledge: L3 Mathematical Structures ? Graph data structure operations
 */

#include "kclique_types.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ================================================================
 * Graph Memory Management
 * ================================================================ */

graph_t *graph_create(int32_t n) {
    if (n <= 0 || n > GRAPH_MAX_VERTICES) {
        fprintf(stderr, "graph_create: invalid n=%d\n", n);
        return NULL;
    }
    graph_t *g = (graph_t *)calloc(1, sizeof(graph_t));
    if (!g) return NULL;

    int64_t matrix_size = (int64_t)n * n;
    g->adj = (uint8_t *)calloc((size_t)matrix_size, sizeof(uint8_t));
    if (!g->adj) {
        free(g);
        return NULL;
    }
    g->n = n;
    g->m = 0;
    g->directed = false;
    g->simple = true;
    return g;
}

void graph_destroy(graph_t *g) {
    if (!g) return;
    free(g->adj);
    g->adj = NULL;
    g->n = 0;
    g->m = 0;
    free(g);
}

graph_t *graph_clone(const graph_t *g) {
    if (!g) return NULL;
    graph_t *copy = graph_create(g->n);
    if (!copy) return NULL;

    int64_t matrix_size = (int64_t)g->n * g->n;
    memcpy(copy->adj, g->adj, (size_t)matrix_size * sizeof(uint8_t));
    copy->m = g->m;
    copy->directed = g->directed;
    copy->simple = g->simple;
    return copy;
}

/* ================================================================
 * Edge Operations
 * ================================================================ */

bool graph_add_edge(graph_t *g, int32_t u, int32_t v) {
    if (!g || u < 0 || u >= g->n || v < 0 || v >= g->n) return false;
    if (g->simple && u == v) return false; /* no self-loops */
    int32_t idx1 = u * g->n + v;
    int32_t idx2 = v * g->n + u;

    if (g->adj[idx1] == 0) {
        g->adj[idx1] = 1;
        if (!g->directed && u != v) {
            g->adj[idx2] = 1;
        }
        g->m++;
        return true;
    }
    return false; /* edge already exists */
}

bool graph_remove_edge(graph_t *g, int32_t u, int32_t v) {
    if (!g || u < 0 || u >= g->n || v < 0 || v >= g->n) return false;
    int32_t idx1 = u * g->n + v;
    int32_t idx2 = v * g->n + u;

    if (g->adj[idx1] == 1) {
        g->adj[idx1] = 0;
        if (!g->directed && u != v) {
            g->adj[idx2] = 0;
        }
        g->m--;
        return true;
    }
    return false;
}

bool graph_has_edge(const graph_t *g, int32_t u, int32_t v) {
    if (!g || u < 0 || u >= g->n || v < 0 || v >= g->n) return false;
    return g->adj[u * g->n + v] != 0;
}

int32_t graph_degree(const graph_t *g, int32_t v) {
    if (!g || v < 0 || v >= g->n) return -1;
    int32_t deg = 0;
    for (int32_t i = 0; i < g->n; i++) {
        if (g->simple && i == v) continue;
        if (g->adj[v * g->n + i]) deg++;
    }
    return deg;
}

int32_t *graph_neighbors(const graph_t *g, int32_t v, int32_t *count) {
    if (!g || !count || v < 0 || v >= g->n) return NULL;
    /* First pass: count neighbors */
    int32_t deg = graph_degree(g, v);
    if (deg < 0) return NULL;
    *count = deg;
    if (deg == 0) return NULL;

    int32_t *neighbors = (int32_t *)malloc((size_t)deg * sizeof(int32_t));
    if (!neighbors) return NULL;

    int32_t idx = 0;
    for (int32_t i = 0; i < g->n; i++) {
        if (g->simple && i == v) continue;
        if (g->adj[v * g->n + i]) {
            neighbors[idx++] = i;
        }
    }
    return neighbors;
}

/* ================================================================
 * Graph Transformations
 * ================================================================ */

graph_t *graph_complement(const graph_t *g) {
    if (!g) return NULL;
    graph_t *comp = graph_create(g->n);
    if (!comp) return NULL;

    for (int32_t i = 0; i < g->n; i++) {
        for (int32_t j = i + 1; j < g->n; j++) {
            if (g->adj[i * g->n + j] == 0) {
                graph_add_edge(comp, i, j);
            }
        }
    }
    comp->directed = g->directed;
    comp->simple = g->simple;
    return comp;
}

graph_t *graph_induced_subgraph(const graph_t *g,
                                 const int32_t *vertices, int32_t count) {
    if (!g || !vertices || count <= 0 || count > g->n) return NULL;
    graph_t *sub = graph_create(count);
    if (!sub) return NULL;

    for (int32_t i = 0; i < count; i++) {
        for (int32_t j = i + 1; j < count; j++) {
            if (g->adj[vertices[i] * g->n + vertices[j]]) {
                graph_add_edge(sub, i, j);
            }
        }
    }
    return sub;
}

/* ================================================================
 * Clique Operations
 * ================================================================ */

void clique_free(clique_t *c) {
    if (!c) return;
    free(c->vertices);
    c->vertices = NULL;
    c->k = 0;
    c->found = false;
}

bool is_adjacent_to_all(const graph_t *g, int32_t v,
                        const int32_t *vertex_set, int32_t set_size) {
    if (!g || !vertex_set || set_size < 0) return false;
    if (set_size == 0) return true;
    for (int32_t i = 0; i < set_size; i++) {
        if (!graph_has_edge(g, v, vertex_set[i])) {
            return false;
        }
    }
    return true;
}

int32_t *common_neighborhood(const graph_t *g,
                              const int32_t *vertex_set, int32_t set_size,
                              int32_t *count) {
    if (!g || !vertex_set || set_size < 0 || !count) {
        if (count) *count = 0;
        return NULL;
    }
    if (set_size == 0) {
        /* Common neighborhood of empty set is all vertices */
        *count = g->n;
        int32_t *all = (int32_t *)malloc((size_t)g->n * sizeof(int32_t));
        if (all) {
            for (int32_t i = 0; i < g->n; i++) all[i] = i;
        }
        return all;
    }

    /* Start with neighbors of first vertex, then intersect */
    int32_t first_count;
    int32_t *common = graph_neighbors(g, vertex_set[0], &first_count);
    if (!common && first_count > 0) return NULL;

    for (int32_t s = 1; s < set_size; s++) {
        int32_t new_count = 0;
        int32_t *temp = (int32_t *)malloc((size_t)first_count * sizeof(int32_t));
        if (!temp) { free(common); return NULL; }
        for (int32_t i = 0; i < first_count; i++) {
            if (graph_has_edge(g, common[i], vertex_set[s])) {
                temp[new_count++] = common[i];
            }
        }
        free(common);
        common = temp;
        first_count = new_count;
        if (first_count == 0) { free(common); *count = 0; return NULL; }
    }
    *count = first_count;
    return common;
}

/* ================================================================
 * Configuration Management
 * ================================================================ */

/**
 * Initialize ETH parameters with default SETH-compliant values.
 * s_3 = 1.0 is the weakest ETH assumption that still implies
 * meaningful lower bounds. Larger s_3 gives stronger bounds.
 *
 * Theorem: If s_3 > 0, then k-Clique not in n^{o(k)} (ETH).
 * If s_k -> 1, then k-Clique not in n^{k-eps} (SETH).
 */
void eth_params_default(eth_params_t *eth) {
    if (!eth) return;
    eth->s_k = ETH_DEFAULT_SK;
    eth->lower_bound = 0.0;
    eth->clique_parameter = 3;
    eth->eth_assumed = true;
}

void seth_params_default(seth_params_t *seth) {
    if (!seth) return;
    seth->cn = SETH_DEFAULT_CN;
    seth->epsilon = 0.01;    /* arbitrary small epsilon */
    seth->sat_k = 3;
    seth->clique_k = 3;
    seth->implied_lower = 3.0 - seth->epsilon;
    seth->seth_assumed = true;
}

conditional_lower_bound_t make_conditional_lower_bound(
    const char *conjecture_name,
    const char *target_problem,
    double exponent,
    bool is_exponential,
    bool is_tight) {
    conditional_lower_bound_t clb;
    clb.conjecture = conjecture_name ? strdup(conjecture_name) : NULL;
    clb.problem = target_problem ? strdup(target_problem) : NULL;
    clb.time_exponent = exponent;
    clb.exponential = is_exponential;
    clb.conditional_constant = 1.0;
    clb.tight = is_tight;
    return clb;
}

void free_clb(conditional_lower_bound_t *clb) {
    if (!clb) return;
    free(clb->conjecture);
    free(clb->problem);
    clb->conjecture = NULL;
    clb->problem = NULL;
}

void print_clb(const conditional_lower_bound_t *clb) {
    if (!clb) { printf("CLB: (null)\n"); return; }
    printf("Conditional Lower Bound:\n");
    printf("  Conjecture: %s\n", clb->conjecture ? clb->conjecture : "(none)");
    printf("  Problem:    %s\n", clb->problem ? clb->problem : "(none)");
    if (clb->exponential) {
        printf("  Bound:      2^{Omega(n)}\n");
    } else {
        printf("  Bound:      n^{%.4f}\n", clb->time_exponent);
    }
    printf("  Tight:      %s\n", clb->tight ? "yes" : "no");
}

double compute_lower_bound_value(int32_t n, double alpha) {
    if (n <= 0 || alpha <= 0.0) return 0.0;
    double val = pow((double)n, alpha);
    if (val > 1e300) return -1.0; /* overflow */
    return val;
}

int compare_clb_hardness(const conditional_lower_bound_t *a,
                         const conditional_lower_bound_t *b) {
    if (!a || !b) return 0;
    if (a->exponential && !b->exponential) return 1;
    if (!a->exponential && b->exponential) return -1;
    if (a->time_exponent > b->time_exponent) return 1;
    if (a->time_exponent < b->time_exponent) return -1;
    return 0;
}

/**
 * Compute binomial coefficient C(n, k) = n! / (k! * (n-k)!).
 * Uses product formula to avoid overflow: C(n,k) = prod_{i=1}^k (n - k + i) / i.
 * Returns -1 on overflow.
 */
int64_t binomial_coefficient(int32_t n, int32_t k) {
    if (k < 0 || k > n) return 0;
    if (k == 0 || k == n) return 1;
    if (k > n - k) k = n - k; /* C(n,k) = C(n,n-k), use smaller k */

    int64_t result = 1;
    for (int32_t i = 1; i <= k; i++) {
        /* Check overflow before multiplication */
        if (result > INT64_MAX / (n - k + i)) return -1;
        result = result * (n - k + i) / i;
    }
    return result;
}

/* next_combination defined in kclique_core.c */
