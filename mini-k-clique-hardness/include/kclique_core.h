/**
 * kclique_core.h ? Core k-Clique operations
 *
 * Fundamental operations for k-Clique checking, enumeration,
 * and maximal clique detection.
 *
 * Reference: Garey & Johnson, "Computers and Intractability" (1979)
 * Reference: Bron & Kerbosch, "Finding all cliques" (CACM 1973)
 * Reference: Tomita, Tanaka, Takahashi, "Worst-case time complexity
 *             for generating all maximal cliques" (TCS 2006)
 *
 * Knowledge: L1 Definitions ? k-Clique problem formalization
 * Knowledge: L3 Mathematical Structures ? Graph operations
 * Knowledge: L5 Algorithms ? Clique detection methods
 * Knowledge: L6 Canonical Problems ? k-Clique, Max-Clique
 */

#ifndef KCLIQUE_CORE_H
#define KCLIQUE_CORE_H

#include "kclique_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ================================================================
 * Core Graph Operations
 * ================================================================ */

/**
 * Allocate and initialize an n-vertex undirected simple graph.
 * All edges initially absent (empty graph on n vertices).
 * Complexity: O(n^2) time, O(n^2) space.
 *
 * @param n  number of vertices (1 <= n <= GRAPH_MAX_VERTICES)
 * @return   pointer to initialized graph_t, or NULL on failure
 */
graph_t *graph_create(int32_t n);

/**
 * Free all resources associated with a graph.
 *
 * @param g  graph to destroy (may be NULL)
 */
void graph_destroy(graph_t *g);

/**
 * Add an undirected edge between vertices u and v.
 * For simple graphs, self-loops are rejected and duplicates are ignored.
 * Complexity: O(1).
 *
 * @param g  graph to modify
 * @param u  first vertex (0 <= u < g->n)
 * @param v  second vertex (0 <= v < g->n)
 * @return   true if edge was added or already existed, false on error
 */
bool graph_add_edge(graph_t *g, int32_t u, int32_t v);

/**
 * Remove an undirected edge between vertices u and v.
 * Complexity: O(1).
 *
 * @param g  graph to modify
 * @param u  first vertex
 * @param v  second vertex
 * @return   true if edge existed and was removed, false otherwise
 */
bool graph_remove_edge(graph_t *g, int32_t u, int32_t v);

/**
 * Check if edge (u,v) exists in the graph.
 * Complexity: O(1).
 *
 * @param g  graph
 * @param u  first vertex
 * @param v  second vertex
 * @return   true if edge exists
 */
bool graph_has_edge(const graph_t *g, int32_t u, int32_t v);

/**
 * Get the degree of vertex v (number of incident edges).
 * Complexity: O(n) for adjacency matrix.
 *
 * @param g  graph
 * @param v  vertex index
 * @return   degree of v, or -1 on error
 */
int32_t graph_degree(const graph_t *g, int32_t v);

/**
 * Get the set of neighbors of vertex v.
 * Caller must free the returned array.
 * Complexity: O(n).
 *
 * @param g      graph
 * @param v      vertex index
 * @param[out] count  number of neighbors returned
 * @return       array of neighbor vertex indices, or NULL on error
 */
int32_t *graph_neighbors(const graph_t *g, int32_t v, int32_t *count);

/**
 * Compute the complement graph: edge (u,v) exists in complement
 * iff it does NOT exist in the original graph (excluding self-loops).
 * This is critical for the k-Clique <-> k-Independent Set duality.
 *
 * @param g  original graph (not modified)
 * @return   new graph = complement(g), or NULL on failure
 */
graph_t *graph_complement(const graph_t *g);

/**
 * Extract the subgraph induced by a set of vertices.
 * Vertices are renumbered 0..count-1 in the induced subgraph.
 *
 * @param g         original graph
 * @param vertices  array of vertex indices to include
 * @param count     number of vertices in the induced subgraph
 * @return          induced subgraph, or NULL on failure
 */
graph_t *graph_induced_subgraph(const graph_t *g,
                                 const int32_t *vertices, int32_t count);

/**
 * Clone a graph (deep copy).
 *
 * @param g  graph to clone
 * @return   independent copy of g, or NULL on failure
 */
graph_t *graph_clone(const graph_t *g);

/* ================================================================
 * k-Clique Decision and Search
 * ================================================================ */

/**
 * Verify that a given set of vertices forms a k-clique in graph g.
 * Checks all k*(k-1)/2 pairs for adjacency.
 * Complexity: O(k^2).
 *
 * @param g         graph
 * @param vertices  array of k vertex indices to check
 * @param k         size of the candidate clique
 * @return          true if vertices form a k-clique
 */
bool is_k_clique(const graph_t *g, const int32_t *vertices, int32_t k);

/**
 * Verify that a given set of vertices forms a clique (any size).
 * Checks pairwise adjacency among all provided vertices.
 * Complexity: O(count^2).
 *
 * @param g         graph
 * @param vertices  array of vertex indices
 * @param count     number of vertices
 * @return          true if the vertices are pairwise adjacent
 */
bool is_clique(const graph_t *g, const int32_t *vertices, int32_t count);

/**
 * Naive brute-force k-clique search: enumerate all C(n,k) subsets
 * and check each for being a clique.
 * Complexity: O(C(n,k) * k^2) = O(n^k) for fixed k.
 * This is the trivial algorithm that SETH suggests is optimal.
 *
 * @param g     graph to search
 * @param k     clique size to find
 * @param[out] result  found clique (caller must free result.vertices if found)
 * @return      true if a k-clique exists
 */
bool find_k_clique_brute(const graph_t *g, int32_t k, clique_t *result);

/**
 * Count the number of k-cliques in the graph.
 * Uses the same brute-force enumeration but counts all occurrences.
 * Complexity: O(C(n,k) * k^2).
 * Note: k-Clique counting is #W[1]-complete (Flum & Grohe 2004).
 *
 * @param g  graph
 * @param k  clique size
 * @return   number of k-cliques, or -1 on error
 */
int64_t count_k_cliques(const graph_t *g, int32_t k);

/**
 * Check if vertex v is adjacent to all vertices in a given set.
 * Complexity: O(|set|).
 *
 * @param g          graph
 * @param v          vertex to check
 * @param vertex_set set of vertices
 * @param set_size   size of the set
 * @return           true if v is adjacent to all vertices in the set
 */
bool is_adjacent_to_all(const graph_t *g, int32_t v,
                        const int32_t *vertex_set, int32_t set_size);

/**
 * Compute the common neighborhood of a set of vertices:
 * { u | forall v in vertex_set: (u,v) in E }.
 * These are vertices that, together with vertex_set, could form
 * a larger clique.
 *
 * @param g          graph
 * @param vertex_set input set of vertices
 * @param set_size   size of the vertex set
 * @param[out] count  number of vertices in common neighborhood
 * @return           array of common neighbor indices (caller must free)
 */
int32_t *common_neighborhood(const graph_t *g,
                              const int32_t *vertex_set, int32_t set_size,
                              int32_t *count);

/* ================================================================
 * Maximum Clique via Bron-Kerbosch
 *
 * The Bron-Kerbosch algorithm (1973) with pivot optimization
 * enumerates all maximal cliques. Modified to also find the
 * maximum clique size.
 *
 * Complexity: O(3^{n/3}) worst-case (Tomita et al. 2006),
 * which is optimal for clique enumeration (Moon-Moser 1965).
 * ================================================================ */

/**
 * Find a maximum clique (largest cardinality) using Bron-Kerbosch
 * with the Tomita pivot heuristic.
 *
 * @param g        graph
 * @param[out] result  maximum clique found (caller must free result.vertices)
 * @return         size k of the maximum clique, or 0 if graph is empty
 */
int32_t max_clique_bron_kerbosch(const graph_t *g, clique_t *result);

/**
 * Enumerate all maximal cliques (inclusion-maximal) in the graph.
 *
 * A maximal clique is a clique that cannot be extended by adding
 * any other vertex.
 *
 * @param g           graph
 * @param[out] cliques  pointer to array of clique_t (caller must free)
 * @param[out] num_cliques  number of maximal cliques found
 * @return            total number enumerated, or -1 on error
 */
int32_t enumerate_maximal_cliques(const graph_t *g,
                                   clique_t **cliques,
                                   int32_t *num_cliques);

/**
 * Free a clique_t structure (frees the internal vertices array).
 */
void clique_free(clique_t *c);

#ifdef __cplusplus
}
#endif

#endif /* KCLIQUE_CORE_H */
