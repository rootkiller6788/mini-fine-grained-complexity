/**
 * kclique_graph.h ? Graph generation and property analysis
 *
 * Functions for generating structured and random graphs,
 * computing graph-theoretic properties relevant to clique analysis.
 *
 * Reference: Erdos & Renyi, "On Random Graphs I" (1959)
 * Reference: Bollobas, "Random Graphs" (2001)
 * Reference: Alon & Spencer, "The Probabilistic Method" (2016)
 * Reference: Chung & Lu, "Complex Graphs and Networks" (2006)
 *
 * Knowledge: L3 Mathematical Structures ? Graph properties
 * Knowledge: L5 Algorithms ? Graph generation, degeneracy
 * Knowledge: L7 Applications ? Random graph models
 */

#ifndef KCLIQUE_GRAPH_H
#define KCLIQUE_GRAPH_H

#include "kclique_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ================================================================
 * Graph Generation
 * ================================================================ */

/**
 * Generate an Erdos-Renyi random graph G(n, p).
 * Each of the n*(n-1)/2 possible edges is included independently
 * with probability p.
 *
 * Properties:
 * - Expected number of edges: p * n*(n-1)/2
 * - Phase transition at p = 1/n for giant component
 * - Expected clique number ~ 2 log_{1/p}(n) for p constant
 * - For p = 1/2, clique number ~ 2 log_2(n)
 *
 * Complexity: O(n^2).
 *
 * @param params  generation parameters
 * @return        generated graph, or NULL on failure
 */
graph_t *graph_generate_erdos_renyi(const er_graph_params_t *params);

/**
 * Generate a complete graph K_n.
 * All n*(n-1)/2 edges are present.
 * Clique number = n.
 *
 * @param n  number of vertices
 * @return   K_n graph
 */
graph_t *graph_generate_complete(int32_t n);

/**
 * Generate an empty graph (independent set) on n vertices.
 * Clique number = 1.
 *
 * @param n  number of vertices
 * @return   empty graph
 */
graph_t *graph_generate_empty(int32_t n);

/**
 * Generate a random k-partite graph with planted k-clique.
 * Useful for testing clique detection algorithms with known ground truth.
 *
 * Partition vertices into k color classes. Within each class, edges
 * exist with probability p_intra. Between classes, edges exist with
 * probability p_inter. Plant one vertex from each class as a clique.
 *
 * @param n            total vertices (must be divisible by k)
 * @param k            number of partitions (and planted clique size)
 * @param p_intra      edge probability within partitions
 * @param p_inter      edge probability between partitions
 * @param seed         random seed
 * @return             generated graph
 */
graph_t *graph_generate_planted_clique(int32_t n, int32_t k,
                                        double p_intra, double p_inter,
                                        int64_t seed);

/**
 * Generate a circulant graph: vertices {0,...,n-1}, edge (i,j)
 * exists iff |i-j| mod n is in the jump set.
 *
 * Circulant graphs have known clique numbers for certain parameter sets.
 *
 * @param n          number of vertices
 * @param jumps      array of jump distances
 * @param num_jumps  number of jump distances
 * @return           generated circulant graph
 */
graph_t *graph_generate_circulant(int32_t n,
                                   const int32_t *jumps,
                                   int32_t num_jumps);

/**
 * Generate the complement of a given graph.
 * Alias for graph_complement (defined in kclique_core.h).
 */
graph_t *graph_generate_complement(const graph_t *g);

/* ================================================================
 * Graph Properties
 * ================================================================ */

/**
 * Compute the degeneracy (d-core number) of a graph.
 * A graph is d-degenerate if every subgraph has a vertex of degree <= d.
 * The degeneracy is the minimum d such that the graph is d-degenerate.
 *
 * Algorithm: repeatedly remove minimum-degree vertex, track max degree seen.
 * Complexity: O(n^2) with adjacency matrix.
 *
 * Significance: any graph with degeneracy d has a (d+1)-coloring,
 * providing an upper bound on clique size: omega <= d + 1.
 *
 * @param g     graph
 * @param[out] ordering  degeneracy ordering of vertices (caller must free)
 * @return      degeneracy value
 */
int32_t graph_degeneracy(const graph_t *g, int32_t **ordering);

/**
 * Compute the core decomposition: for each k, the k-core is the
 * maximal subgraph where all vertices have degree >= k.
 *
 * Stores core values: core_value[v] = max k such that v is in k-core.
 *
 * @param g           graph
 * @param[out] core_values  core number for each vertex (caller allocates n ints)
 */
void graph_core_decomposition(const graph_t *g, int32_t *core_values);

/**
 * Compute the degree distribution of the graph.
 *
 * @param g              graph
 * @param[out] deg_count  array of size g->n+1 where deg_count[d] counts
 *                        vertices of degree d
 * @param[out] max_deg    maximum degree observed
 */
void graph_degree_distribution(const graph_t *g,
                                int32_t *deg_count,
                                int32_t *max_deg);

/**
 * Compute the clustering coefficient of vertex v:
 * (number of edges between neighbors of v) / (d(v) choose 2).
 *
 * @param g  graph
 * @param v  vertex index
 * @return   clustering coefficient in [0,1], or -1 on error
 */
double graph_clustering_coefficient(const graph_t *g, int32_t v);

/**
 * Compute the average clustering coefficient of the graph.
 *
 * @param g  graph
 * @return   average clustering coefficient
 */
double graph_average_clustering(const graph_t *g);

/**
 * Check if the graph is triangle-free (contains no 3-clique).
 * This is a special case of k-Clique for k=3.
 * Complexity: O(n*d_max) naive, O(n^omega) via matrix multiplication.
 *
 * @param g  graph
 * @return   true if the graph has no triangles
 */
bool graph_is_triangle_free(const graph_t *g);

/**
 * Compute the Turan bound: maximum number of edges in a K_{r+1}-free
 * graph on n vertices is (1 - 1/r) * n^2 / 2 (Turan's theorem, 1941).
 * This gives a necessary condition for clique existence.
 *
 * @param n      number of vertices
 * @param r      forbidden clique size minus 1
 * @return       Turan bound on edge count
 */
int64_t turan_bound_edges(int32_t n, int32_t r);

/**
 * Check if the graph exceeds the Turan bound for K_{k}-free graphs.
 * If |E| > turan_bound(n, k-1), then the graph MUST contain a k-clique
 * (by Turan's theorem).
 *
 * @param g  graph
 * @param k  clique size to test
 * @return   true if a k-clique is guaranteed to exist (by Turan)
 */
bool turan_guarantees_clique(const graph_t *g, int32_t k);

/**
 * Compute the density of the graph: 2|E| / (n(n-1)).
 *
 * @param g  graph
 * @return   density in [0, 1]
 */
double graph_density(const graph_t *g);

#ifdef __cplusplus
}
#endif

#endif /* KCLIQUE_GRAPH_H */
