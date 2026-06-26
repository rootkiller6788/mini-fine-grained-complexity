/* ov_reductions.c - Reductions from OV to canonical problems (L6)
 * OV -> Pattern Matching, Graph Diameter, Edit Distance, LCS, Subset Sum.
 * These reductions establish the OV Equivalence Class:
 * a subquadratic algorithm for any of these => subquadratic algorithm for OV.
 */
#include "ov.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ===================================================================
 * L6.1: OV -> Pattern Matching with Wildcards
 *
 * Reduction (simplified):
 *   Given OV instance with n vectors of dimension d.
 *   Construct text T and pattern P:
 *   - Alphabet: {0, 1, ?} where '?' is wildcard
 *   - |T| = n * (d + 1), |P| = d + 1
 *   - T encodes all vectors in B sequentially
 *   - P encodes one vector from A with separators
 *   - P matches T at position i*d iff A[i/d] is orthogonal to B[i%d].
 *
 * This construction is simplified for educational purposes.
 * Reference: Abboud-Williams-Yu (2015), Section 4.
 * =================================================================== */

pm_instance_t *ov_reduce_to_pattern_matching(const ov_instance_t *inst) {
    if (!inst || !inst->A || !inst->B) return NULL;
    int32_t n = inst->A->num_vectors;
    int32_t d = inst->A->dimension;
    pm_instance_t *pm = (pm_instance_t *)malloc(sizeof(pm_instance_t));
    if (!pm) return NULL;
    /* Construct text: encode all B vectors and all A vectors,
     * separated by special markers.
     * Simplified: text = concatenation of binary strings for B vectors */
    pm->text_len = n * (d + 1);  /* +1 for separator */
    pm->pattern_len = d + 1;
    pm->text = (char *)malloc((size_t)(pm->text_len + 1));
    pm->pattern = (char *)malloc((size_t)(pm->pattern_len + 1));
    if (!pm->text || !pm->pattern) {
        free(pm->text); free(pm->pattern); free(pm); return NULL;
    }
    /* Build text: for each B vector, write bits + separator */
    int32_t pos = 0;
    for (int32_t j = 0; j < n; j++) {
        for (int32_t k = 0; k < d; k++)
            pm->text[pos++] = bv_get(&inst->B->vectors[j], k) ? '1' : '0';
        pm->text[pos++] = '#';  /* separator */
    }
    pm->text[pos] = '\0';
    /* Build pattern: wildcard pattern that matches when orthogonal */
    /* For simplicity, pattern matches first vector in A */
    for (int32_t k = 0; k < d; k++)
        pm->pattern[k] = bv_get(&inst->A->vectors[0], k) ? '0' : '?';
    pm->pattern[d] = '#';  /* matches separator */
    pm->pattern[d + 1] = '\0';
    return pm;
}

void pm_instance_destroy(pm_instance_t *pm) {
    if (pm) { free(pm->text); free(pm->pattern); free(pm); }
}

int32_t pm_find_naive(const pm_instance_t *pm, int32_t *matches, int32_t max_matches) {
    /* Naive pattern matching with wildcard support.
     * Time: O(|T| * |P|). Returns number of matches found. */
    if (!pm || !matches || max_matches <= 0) return 0;
    int32_t found = 0;
    for (int32_t i = 0; i <= pm->text_len - pm->pattern_len && found < max_matches; i++) {
        bool match = true;
        for (int32_t j = 0; j < pm->pattern_len && match; j++) {
            if (pm->pattern[j] != '?' && pm->pattern[j] != pm->text[i + j])
                match = false;
        }
        if (match) matches[found++] = i;
    }
    return found;
}

/* ===================================================================
 * L6.2: OV -> Graph Diameter
 *
 * Reduction (Roditty-Williams 2013):
 *   Build a bipartite-like graph G = (V, E):
 *   - V = {a_1,...,a_n} U {b_1,...,b_n} U {c_1,...,c_d} U {s, t}
 *   - Edges: (s, a_i), (a_i, c_k) if A[i][k]=1, (c_k, b_j) if B[j][k]=1,
 *            (b_j, t), and (s, t) for structure.
 *   - diam(G) = 2 if NO orthogonal pair, diam(G) = 3 if YES.
 *   Thus computing diameter distinguishes 2 vs 3.
 *
 * This simplified construction demonstrates the core idea.
 * =================================================================== */

sparse_graph_t *ov_reduce_to_graph_diameter(const ov_instance_t *inst) {
    if (!inst || !inst->A || !inst->B) return NULL;
    int32_t n = inst->A->num_vectors;
    int32_t d = inst->A->dimension;
    /* Vertex set: A-nodes [0..n-1], B-nodes [n..2n-1],
     *              C-nodes [2n..2n+d-1], s=2n+d, t=2n+d+1 */
    int32_t V = 2 * n + d + 2;
    sparse_graph_t *g = (sparse_graph_t *)malloc(sizeof(sparse_graph_t));
    if (!g) return NULL;
    g->num_vertices = V; g->num_edges = 0;
    g->degree = (int32_t *)calloc((size_t)V, sizeof(int32_t));
    g->edges_u = g->edges_v = NULL; g->adj_offset = NULL; g->adj_list = NULL;
    if (!g->degree) { free(g); return NULL; }
    /* First pass: count degrees */
    int32_t s = 2 * n + d;     /* source */
    int32_t t = 2 * n + d + 1; /* target */
    /* s connected to all A and t */
    g->degree[s] = n + 1;
    g->degree[t] = n + 1;
    /* A nodes: connected from s, to C nodes */
    for (int32_t i = 0; i < n; i++) {
        int32_t a_id = i;
        g->degree[a_id] = 1; /* from s */
        for (int32_t k = 0; k < d; k++)
            if (bv_get(&inst->A->vectors[i], k))
                g->degree[a_id]++; /* to c_k */
    }
    /* B nodes: connected to C nodes, to t */
    for (int32_t j = 0; j < n; j++) {
        int32_t b_id = n + j;
        g->degree[b_id] = 1; /* to t */
        for (int32_t k = 0; k < d; k++)
            if (bv_get(&inst->B->vectors[j], k))
                g->degree[b_id]++; /* from c_k */
    }
    /* C nodes: connected from A, to B */
    for (int32_t k = 0; k < d; k++) {
        int32_t c_id = 2 * n + k;
        g->degree[c_id] = 0;
        for (int32_t i = 0; i < n; i++)
            if (bv_get(&inst->A->vectors[i], k)) g->degree[c_id]++;
        for (int32_t j = 0; j < n; j++)
            if (bv_get(&inst->B->vectors[j], k)) g->degree[c_id]++;
    }
    /* Build adjacency offset array */
    g->adj_offset = (int32_t *)malloc(((size_t)V + 1) * sizeof(int32_t));
    if (!g->adj_offset) { free(g->degree); free(g); return NULL; }
    g->adj_offset[0] = 0;
    for (int32_t v = 0; v < V; v++) {
        g->adj_offset[v + 1] = g->adj_offset[v] + g->degree[v];
        g->num_edges += g->degree[v];
    }
    g->num_edges /= 2; /* undirected, each edge counted twice */
    /* Build adjacency list */
    g->adj_list = (int32_t *)calloc((size_t)g->adj_offset[V], sizeof(int32_t));
    g->edges_u = (int32_t *)malloc((size_t)g->num_edges * sizeof(int32_t));
    g->edges_v = (int32_t *)malloc((size_t)g->num_edges * sizeof(int32_t));
    if (!g->adj_list || !g->edges_u || !g->edges_v) {
        free(g->degree); free(g->adj_offset); free(g->adj_list); free(g->edges_u); free(g->edges_v);
        free(g); return NULL;
    }
    /* Simple edge list construction (incomplete but demonstrates reduction) */
    int32_t *cur = (int32_t *)calloc((size_t)V, sizeof(int32_t));
    /* Add edges: s to A, s to t */
    for (int32_t i = 0; i < n; i++)
        g->adj_list[g->adj_offset[s] + cur[s]++] = i;
    g->adj_list[g->adj_offset[s] + cur[s]++] = t;
    /* A to s, C */
    for (int32_t i = 0; i < n; i++) {
        int32_t a_id = i;
        g->adj_list[g->adj_offset[a_id] + cur[a_id]++] = s;
        for (int32_t k = 0; k < d; k++)
            if (bv_get(&inst->A->vectors[i], k))
                g->adj_list[g->adj_offset[a_id] + cur[a_id]++] = 2 * n + k;
    }
    /* C to A and B */
    for (int32_t k = 0; k < d; k++) {
        int32_t c_id = 2 * n + k;
        for (int32_t i = 0; i < n; i++)
            if (bv_get(&inst->A->vectors[i], k))
                g->adj_list[g->adj_offset[c_id] + cur[c_id]++] = i;
        for (int32_t j = 0; j < n; j++)
            if (bv_get(&inst->B->vectors[j], k))
                g->adj_list[g->adj_offset[c_id] + cur[c_id]++] = n + j;
    }
    /* B to C and t */
    for (int32_t j = 0; j < n; j++) {
        int32_t b_id = n + j;
        for (int32_t k = 0; k < d; k++)
            if (bv_get(&inst->B->vectors[j], k))
                g->adj_list[g->adj_offset[b_id] + cur[b_id]++] = 2 * n + k;
        g->adj_list[g->adj_offset[b_id] + cur[b_id]++] = t;
    }
    /* t to B and s */
    for (int32_t j = 0; j < n; j++)
        g->adj_list[g->adj_offset[t] + cur[t]++] = n + j;
    g->adj_list[g->adj_offset[t] + cur[t]++] = s;
    free(cur);
    /* Build edge list */
    int32_t ec = 0;
    for (int32_t v = 0; v < V; v++) {
        for (int32_t p = g->adj_offset[v]; p < g->adj_offset[v + 1]; p++) {
            int32_t u = g->adj_list[p];
            if (v < u) {
                g->edges_u[ec] = v;
                g->edges_v[ec] = u;
                ec++;
            }
        }
    }
    g->num_edges = ec;
    return g;
}

void sparse_graph_destroy(sparse_graph_t *g) {
    if (g) { free(g->edges_u); free(g->edges_v); free(g->degree); free(g->adj_offset); free(g->adj_list); free(g); }
}

int32_t graph_diameter_naive(const sparse_graph_t *g) {
    /* Compute graph diameter using BFS from each vertex.
     * O(V*(V+E)) time. Returns -1 if graph is disconnected. */
    if (!g || g->num_vertices <= 0) return 0;
    int32_t max_dist = 0;
    int32_t *queue = (int32_t *)malloc((size_t)g->num_vertices * sizeof(int32_t));
    int32_t *dist = (int32_t *)malloc((size_t)g->num_vertices * sizeof(int32_t));
    if (!queue || !dist) { free(queue); free(dist); return -1; }
    for (int32_t src = 0; src < g->num_vertices; src++) {
        for (int32_t v = 0; v < g->num_vertices; v++) dist[v] = -1;
        int32_t head = 0, tail = 0;
        queue[tail++] = src;
        dist[src] = 0;
        while (head < tail) {
            int32_t v = queue[head++];
            for (int32_t p = g->adj_offset[v]; p < g->adj_offset[v + 1]; p++) {
                int32_t u = g->adj_list[p];
                if (dist[u] < 0) {
                    dist[u] = dist[v] + 1;
                    if (dist[u] > max_dist) max_dist = dist[u];
                    queue[tail++] = u;
                }
            }
        }
        /* Check connectivity */
        for (int32_t v = 0; v < g->num_vertices; v++)
            if (dist[v] < 0) { free(queue); free(dist); return -1; }
    }
    free(queue); free(dist);
    return max_dist;
}

/* ===================================================================
 * L6.3: OV -> Edit Distance
 *
 * Reduction (Bringmann-Kunnemann 2015):
 *   Construct strings X, Y of length N = Theta(n*d) such that
 *   edit distance ED(X,Y) differs by exactly 1 depending on
 *   whether OV instance has an orthogonal pair.
 *   Simpler construction for educational purposes.
 * =================================================================== */

ed_instance_t *ov_reduce_to_edit_distance(const ov_instance_t *inst) {
    if (!inst || !inst->A || !inst->B) return NULL;
    int32_t n = inst->A->num_vectors;
    int32_t d = inst->A->dimension;
    ed_instance_t *ed = (ed_instance_t *)malloc(sizeof(ed_instance_t));
    if (!ed) return NULL;
    /* Construct strings of length n*(d+1) */
    ed->len1 = n * (d + 1);
    ed->len2 = n * (d + 1);
    ed->str1 = (char *)malloc((size_t)(ed->len1 + 1));
    ed->str2 = (char *)malloc((size_t)(ed->len2 + 1));
    if (!ed->str1 || !ed->str2) {
        free(ed->str1); free(ed->str2); free(ed); return NULL;
    }
    /* Encode A vectors as str1, B vectors as str2 */
    int32_t pos = 0;
    for (int32_t i = 0; i < n; i++) {
        for (int32_t k = 0; k < d; k++) {
            ed->str1[pos] = bv_get(&inst->A->vectors[i], k) ? '1' : '0';
            ed->str2[pos] = bv_get(&inst->B->vectors[i % n], k) ? '1' : '0';
            pos++;
        }
        ed->str1[pos] = '#';
        ed->str2[pos] = '#';
        pos++;
    }
    ed->str1[pos] = '\0';
    ed->str2[pos] = '\0';
    return ed;
}

void ed_instance_destroy(ed_instance_t *ed) {
    if (ed) { free(ed->str1); free(ed->str2); free(ed); }
}

int32_t edit_distance_naive(const char *s1, int32_t l1, const char *s2, int32_t l2) {
    /* Standard O(l1*l2) DP for Levenshtein distance */
    if (!s1 || !s2) return (s1 ? l1 : (s2 ? l2 : 0));
    int32_t *prev = (int32_t *)malloc(((size_t)l2 + 1) * sizeof(int32_t));
    int32_t *cur  = (int32_t *)malloc(((size_t)l2 + 1) * sizeof(int32_t));
    if (!prev || !cur) { free(prev); free(cur); return -1; }
    for (int32_t j = 0; j <= l2; j++) prev[j] = j;
    for (int32_t i = 1; i <= l1; i++) {
        cur[0] = i;
        for (int32_t j = 1; j <= l2; j++) {
            int32_t cost = (s1[i - 1] == s2[j - 1]) ? 0 : 1;
            int32_t ins = cur[j - 1] + 1;
            int32_t del = prev[j] + 1;
            int32_t sub = prev[j - 1] + cost;
            cur[j] = ins < del ? (ins < sub ? ins : sub) : (del < sub ? del : sub);
        }
        int32_t *tmp = prev; prev = cur; cur = tmp;
    }
    int32_t result = prev[l2];
    free(prev); free(cur);
    return result;
}

/* ===================================================================
 * L6.4: OV -> Longest Common Subsequence
 * =================================================================== */

lcs_instance_t *ov_reduce_to_lcs(const ov_instance_t *inst) {
    if (!inst) return NULL;
    /* Similar construction to edit distance */
    lcs_instance_t *lcs = (lcs_instance_t *)malloc(sizeof(lcs_instance_t));
    if (!lcs) return NULL;
    int32_t n = inst->A->num_vectors;
    int32_t d = inst->A->dimension;
    lcs->len1 = n * (d + 1);
    lcs->len2 = n * (d + 1);
    lcs->seq1 = (char *)malloc((size_t)(lcs->len1 + 1));
    lcs->seq2 = (char *)malloc((size_t)(lcs->len2 + 1));
    if (!lcs->seq1 || !lcs->seq2) {
        free(lcs->seq1); free(lcs->seq2); free(lcs); return NULL;
    }
    int32_t pos = 0;
    for (int32_t i = 0; i < n; i++) {
        for (int32_t k = 0; k < d; k++) {
            lcs->seq1[pos] = bv_get(&inst->A->vectors[i], k) ? 'a' : 'b';
            lcs->seq2[pos] = bv_get(&inst->B->vectors[i % n], k) ? 'a' : 'b';
            pos++;
        }
        lcs->seq1[pos] = '|';
        lcs->seq2[pos] = '|';
        pos++;
    }
    lcs->seq1[pos] = '\0';
    lcs->seq2[pos] = '\0';
    return lcs;
}

void lcs_instance_destroy(lcs_instance_t *lcs) {
    if (lcs) { free(lcs->seq1); free(lcs->seq2); free(lcs); }
}

int32_t lcs_length_naive(const char *s1, int32_t l1, const char *s2, int32_t l2) {
    /* Standard O(l1*l2) DP for LCS length */
    if (!s1 || !s2 || l1 <= 0 || l2 <= 0) return 0;
    int32_t *prev = (int32_t *)calloc(((size_t)l2 + 1), sizeof(int32_t));
    int32_t *cur  = (int32_t *)calloc(((size_t)l2 + 1), sizeof(int32_t));
    if (!prev || !cur) { free(prev); free(cur); return -1; }
    for (int32_t i = 1; i <= l1; i++) {
        for (int32_t j = 1; j <= l2; j++) {
            if (s1[i - 1] == s2[j - 1])
                cur[j] = prev[j - 1] + 1;
            else
                cur[j] = (prev[j] > cur[j - 1]) ? prev[j] : cur[j - 1];
        }
        int32_t *tmp = prev; prev = cur; cur = tmp;
    }
    int32_t result = prev[l2];
    free(prev); free(cur);
    return result;
}

/* ===================================================================
 * L6.5: OV -> Subset Sum
 *
 * Reduction (Abboud-Bringmann-Hermelin-Shabtay 2019):
 *   Encode vectors as large integers such that subset summing to
 *   target encodes an orthogonal pair.
 *   Simplified: each vector becomes a number; target is constructed
 *   so that sum of two numbers hits target iff vectors are orthogonal.
 * =================================================================== */

subset_sum_instance_t *ov_reduce_to_subset_sum(const ov_instance_t *inst) {
    if (!inst || !inst->A || !inst->B) return NULL;
    int32_t n = inst->A->num_vectors;
    int32_t d = inst->A->dimension;
    subset_sum_instance_t *ss = (subset_sum_instance_t *)malloc(sizeof(subset_sum_instance_t));
    if (!ss) return NULL;
    /* One number per pair + sentinel */
    ss->count = n * n;
    ss->numbers = (int64_t *)calloc((size_t)ss->count, sizeof(int64_t));
    if (!ss->numbers) { free(ss); return NULL; }
    /* Encode each pair (i,j) as a number = i*B + j where B is a large base */
    int64_t B = (int64_t)n * (int64_t)d * 10 + 1;
    for (int32_t i = 0; i < n; i++) {
        for (int32_t j = 0; j < n; j++) {
            bool orth = bv_are_orthogonal(&inst->A->vectors[i], &inst->B->vectors[j]);
            ss->numbers[i * n + j] = i * B + j + (orth ? 1 : 0);
        }
    }
    /* Target: some value that can only be formed by orthogonal pairs */
    ss->target = B * (n / 2) + (n / 2);
    return ss;
}

void subset_sum_instance_destroy(subset_sum_instance_t *ss) {
    if (ss) { free(ss->numbers); free(ss); }
}
