#ifndef FINE_GRAINED_REDUCTION_H
#define FINE_GRAINED_REDUCTION_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <math.h>
#include "equiv_classes.h"

/*
 * fine_grained_reduction.h -- Fine-Grained Reductions Framework
 *
 * A fine-grained reduction is a reduction that runs in time strictly
 * less than the conjectured lower bound for the target class.
 * This allows us to transfer lower bounds: if problem A is conjectured
 * to require time T(n) and A reduces to B in time o(T(n)), then B
 * also requires time T(n) (conditionally).
 *
 * Types of fine-grained reductions:
 *   1. Subcubic reductions: O(n^{3-epsilon}) for some epsilon > 0.
 *      Used for the subcubic equivalence class.
 *   2. Subquadratic reductions: O(n^{2-epsilon}) for some epsilon > 0.
 *      Used for the subquadratic and 3SUM equivalence classes.
 *   3. Linear / near-linear reductions: O(n polylog n).
 *      Stronger form of fine-grained reduction.
 *
 * Key Properties (L4):
 *   - Transitivity: If A reduces to B and B reduces to C under
 *     fine-grained reductions, then A reduces to C.
 *   - Composition: The composition of two fine-grained reductions
 *     is a fine-grained reduction (may accumulate overhead).
 *
 * References:
 *   Williams (2015), "Hardness of Easy Problems", ICDT
 *   Bringmann (2019), "Fine-Grained Complexity Theory", SIGACT News
 *
 * L1: Fine-grained reduction types
 * L2: Reduction-preserving equivalences
 * L3: Graph transformations, polynomial encoding
 * L4: Transitivity, composition, completeness
 * L5: Specific reduction constructions
 * L6: SETH-based reductions, APSP-based reductions, 3SUM-based reductions
 */

/* ---- L1: Reduction Graph Node ---- */

/*
 * A node in the fine-grained reduction graph.
 * Each node represents a computational problem.
 */
typedef struct fg_node_t {
    problem_id_t        id;
    char                name[80];
    equiv_class_id_t    class_id;
    double              best_exponent;
    double              conjectured_lower_bound;
    int32_t             num_outgoing;       /* Number of outgoing reductions */
    int32_t             num_incoming;       /* Number of incoming reductions */
    struct fg_edge_t   *out_edges;          /* Linked list of outgoing edges */
    struct fg_edge_t   *in_edges;           /* Linked list of incoming edges */
} fg_node_t;

/* ---- L1: Reduction Graph Edge ---- */

/*
 * An edge in the fine-grained reduction graph.
 * Directed from source problem to target problem.
 */
typedef struct fg_edge_t {
    problem_id_t        from_id;
    problem_id_t        to_id;
    fg_reduction_type_t type;
    double              exponent;           /* Reduction runs in O(n^exponent) */
    int32_t             blowup;             /* Output instance size = O(n^blowup) */
    bool                is_tight;           /* Is the reduction optimal? */
    char                citation[160];
    struct fg_edge_t   *next_out;           /* Next outgoing edge from the same source */
    struct fg_edge_t   *next_in;            /* Next incoming edge to the same target */
} fg_edge_t;

/* ---- L1: Fine-Grained Reduction Graph ---- */

typedef struct {
    fg_node_t  **nodes;         /* Array of node pointers */
    int32_t      num_nodes;
    int32_t      capacity;
    fg_edge_t  **edges;         /* Array of edge pointers */
    int32_t      num_edges;
    int32_t      edge_capacity;
} fg_reduction_graph_t;

/* ---- L2: Reduction Graph Operations ---- */

fg_reduction_graph_t *fg_graph_create(int32_t node_cap, int32_t edge_cap);
void fg_graph_free(fg_reduction_graph_t *g);

/*
 * Add a problem node to the reduction graph.
 * Returns the node pointer, or NULL on error.
 */
fg_node_t *fg_graph_add_node(fg_reduction_graph_t *g, const char *name,
                              equiv_class_id_t class_id, double best_exp,
                              double conj_lb);

/*
 * Add a fine-grained reduction edge between two problems.
 * The edge encodes: "from reduces to to" with the given parameters.
 * Returns the edge pointer, or NULL on error.
 */
fg_edge_t *fg_graph_add_edge(fg_reduction_graph_t *g,
                              problem_id_t from, problem_id_t to,
                              fg_reduction_type_t type, double exponent,
                              int32_t blowup, bool is_tight,
                              const char *citation);

/*
 * Find a problem node by name. Returns NULL if not found.
 */
fg_node_t *fg_graph_find_node(const fg_reduction_graph_t *g, const char *name);

/* ---- L2: Reduction Chain Computation ---- */

/*
 * Compute the shortest (minimum total exponent overhead) reduction
 * chain from problem 'from' to problem 'to'.
 *
 * Uses Dijkstra-like algorithm on the reduction graph where edge
 * weights are the exponent overheads.
 *
 * Returns the minimum chain exponent sum, or -1 if no path exists.
 * Fills 'path' and 'path_len' with the sequence of problem IDs
 * along the shortest path (if path is non-NULL).
 */
double fg_shortest_reduction_chain(const fg_reduction_graph_t *g,
                                    problem_id_t from, problem_id_t to,
                                    problem_id_t *path, int32_t *path_len,
                                    int32_t max_path_len);

/*
 * Compute the transitive closure of the reduction graph.
 * For each pair of problems (u,v), determine if v is reachable from u
 * and record the minimum exponent overhead along any path.
 *
 * Results stored in a |V| x |V| matrix: closure[u][v] = min exponent sum.
 * -1.0 means unreachable.
 */
double **fg_transitive_closure(const fg_reduction_graph_t *g);

/*
 * Free the transitive closure matrix.
 */
void fg_transitive_closure_free(double **closure, int32_t n);

/* ---- L2: Equivalence Class from Reduction Graph ---- */

/*
 * Given a reduction graph, identify all equivalence classes.
 * Two problems are in the same class if they are mutually reachable
 * via fine-grained reductions.
 *
 * This computes the strongly connected components (SCCs) of the
 * reduction graph, using Tarjan's algorithm.
 *
 * Returns an array mapping each node ID to its SCC representative.
 * 'num_classes' is set to the number of distinct classes found.
 */
problem_id_t *fg_find_equivalence_classes(const fg_reduction_graph_t *g,
                                           int32_t *num_classes);

/*
 * Check if the reduction graph is "sound" for a given class:
 * all problems in the class are mutually reducible under reductions
 * faster than the class threshold.
 *
 * Returns true if the class forms a clique of mutual fine-grained reductions.
 */
bool fg_class_is_sound(const fg_reduction_graph_t *g, equiv_class_id_t cid);

/* ---- L4: Composition of Fine-Grained Reductions ---- */

/*
 * Compose two fine-grained reductions into one.
 * If P reduces to Q with exponent a and blowup b1, and
 * Q reduces to R with exponent b and blowup b2 (on instances of size
 * n' = O(n^{b1})), then P reduces to R with exponent
 * max(a, b * b1) and total blowup b1 * b2.
 *
 * Returns the composed edge parameters.
 * out_exponent receives the composed runtime exponent.
 * out_blowup receives the composed blowup factor.
 */
bool fg_compose_reductions(double exponent1, int32_t blowup1,
                            double exponent2, int32_t blowup2,
                            double *out_exponent, int32_t *out_blowup);

/*
 * Compose a chain of reductions and compute the total exponent.
 * The chain is given as an array of edge pointers in order.
 *
 * Returns the total exponent for reducing the first problem
 * to the last, or -1 if composition fails.
 */
double fg_compose_chain(fg_edge_t **chain, int32_t chain_len);

/* ---- L5: Specific Reduction Constructions ---- */

/*
 * SETH-based reduction: k-SAT to OV.
 * Theorem (Williams 2005): For any epsilon > 0, if OV can be solved
 * in O(n^{2-epsilon}) time, then SETH is false.
 *
 * This function encodes the reduction as an edge in the graph.
 * Given a k-SAT instance with n vars, produces an OV instance with
 * N = 2^{n/2} vectors of dimension d = O(n).
 */
bool fg_seth_to_ov_reduction(double epsilon, int32_t k,
                              int32_t *out_n, int32_t *out_d,
                              double *out_exponent);

/*
 * OV to Edit Distance reduction.
 * Theorem (Backurs & Indyk 2016): If Edit Distance can be solved in
 * O(n^{2-epsilon}) time, then OV conjecture is false (so SETH false).
 *
 * Encodes two sets of n boolean vectors as two strings of length O(n*d).
 */
bool fg_ov_to_edit_distance_reduction(int32_t n, int32_t d,
                                       int32_t *out_str_len,
                                       double *out_exponent);

/* ---- L7: Applications of Fine-Grained Reductions ---- */

/*
 * Classify a new problem by finding the tightest known lower bounds
 * via reductions from canonical hard problems.
 *
 * Returns the ID of the hardest problem known to reduce to the new problem,
 * along with the reduction exponent.
 */
problem_id_t fg_classify_new_problem(const fg_reduction_graph_t *g,
                                      const char *name, double current_best_exponent,
                                      double *hardness_exponent);

/*
 * Print the fine-grained complexity landscape:
 * all equivalence classes and their relationships.
 */
void fg_print_landscape(const fg_reduction_graph_t *g);

/*
 * Generate a DOT graph description of the reduction graph
 * for visualization with Graphviz.
 */
void fg_print_dot(const fg_reduction_graph_t *g, const char *filename);

/* ---- L8: Advanced Topics ---- */

/*
 * Check if a reduction is "barrier-breaking":
 * would the existence of a sub-threshold reduction from A to B
 * refute a major complexity conjecture?
 *
 * For example, a subcubic reduction from APSP to an O(n^2)-time problem
 * would break the cubic barrier, which is considered highly unlikely.
 */
bool fg_is_barrier_breaking(const fg_reduction_graph_t *g, problem_id_t from,
                             problem_id_t to, double exponent);

/*
 * Compute the "fine-grained diameter" of the entire problem landscape:
 * the maximum distance (in exponent overhead) between any two problems
 * that are in the same equivalence class.
 */
double fg_landscape_diameter(const fg_reduction_graph_t *g);

#endif /* FINE_GRAINED_REDUCTION_H */