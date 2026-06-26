/**
 * kclique_reduction.h ? Reductions involving k-Clique
 *
 * Parameterized and classical reductions:
 * - k-Clique <-> k-Independent Set (via complement)
 * - k-Clique -> k-Vertex Cover (via complement + size)
 * - 3SAT -> k-Clique (FPT reduction, W[1]-hardness proof)
 * - k-Clique -> k x k Clique (matrix/tensor formulation)
 *
 * Reference: Downey & Fellows, "Parameterized Complexity" (1999), Ch. 7
 * Reference: Flum & Grohe, "Parameterized Complexity Theory" (2006)
 * Reference: Chen & Flum, "Parameterized Complexity and Subexponential Time" (2007)
 * Reference: Garey & Johnson, "Computers and Intractability" (1979)
 *
 * Knowledge: L4 Fundamental Laws ? W[1]-completeness of k-Clique
 * Knowledge: L6 Canonical Problems ? k-IS, k-VC, k-Clique reductions
 * Knowledge: L8 Advanced Topics ? FPT reduction properties
 */

#ifndef KCLIQUE_REDUCTION_H
#define KCLIQUE_REDUCTION_H

#include "kclique_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ================================================================
 * k-Clique <-> k-Independent Set Duality
 *
 * An independent set in G is a clique in complement(G).
 * This is a parameter-preserving reduction (k' = k) and
 * runs in polynomial time O(n^2).
 *
 * Since k-Independent Set is W[1]-complete, so is k-Clique.
 * ================================================================ */

/**
 * Convert a k-Clique instance to a k-Independent Set instance.
 * Simply returns the complement graph.
 * Runtime: O(n^2), parameter k preserved.
 *
 * @param g  input graph for k-Clique
 * @return    complement graph for k-Independent Set
 */
graph_t *kclique_to_kindependentset(const graph_t *g);

/**
 * Convert a k-Independent Set instance to a k-Clique instance.
 * Identical to kclique_to_kindependentset (complement is self-inverse).
 *
 * @param g  input graph for k-Independent Set
 * @return   complement graph for k-Clique
 */
graph_t *kindependentset_to_kclique(const graph_t *g);

/* ================================================================
 * k-Clique -> k-Vertex Cover Reduction
 *
 * A graph G has a k-clique iff its complement has an (n-k)-vertex cover.
 * This is a polynomial-time reduction with parameter change k' = n - k.
 * ================================================================ */

/**
 * Reduce k-Clique to (n-k)-Vertex Cover.
 * G has a k-clique iff complement(G) has a vertex cover of size n-k.
 *
 * @param g          input graph
 * @param k          clique size parameter
 * @param[out] vc_k   vertex cover parameter (n - k)
 * @return           complement graph for the vertex cover instance
 */
graph_t *kclique_to_vertexcover(const graph_t *g, int32_t k, int32_t *vc_k);

/* ================================================================
 * 3SAT -> k-Clique FPT Reduction (W[1]-hardness)
 *
 * This is the canonical FPT reduction proving k-Clique is W[1]-hard.
 * Given a 3SAT formula with m clauses, construct a graph with
 * m "column" groups, each containing vertices corresponding to
 * satisfying assignments of clauses.
 *
 * Construction:
 * - For each clause C_j, create 7 vertices (one per satisfying
 *   partial assignment to the 3 literals in C_j)
 * - Connect two vertices u, v iff:
 *   1. They come from different clauses
 *   2. Their partial assignments are compatible (no contradictory literals)
 * - The formula is satisfiable iff the graph has an m-clique
 *
 * This is an FPT reduction: runtime is polynomial in |phi|,
 * and the parameter k' = m is the number of clauses.
 * ================================================================ */

/**
 * Datum: one variable assignment in a clause.
 * variable_index: 0-based variable index
 * value: true for positive literal, false for negative
 */
typedef struct {
    int32_t variable_index;
    bool    value;
} literal_assignment_t;

/**
 * Represent a 3SAT instance for reduction to k-Clique.
 */
typedef struct {
    int32_t  num_variables;       /* number of Boolean variables */
    int32_t  num_clauses;         /* number of clauses (each has 3 literals) */
    int32_t *clause_literals;     /* 3*num_clauses entries:
                                     clause_literals[3*c + i] = variable index,
                                     sign encoded via bit 31:
                                     positive if (literal & 0x80000000) == 0,
                                     negative if (literal & 0x80000000) != 0 */
} sat3_instance_t;

/**
 * Construct the FPT reduction graph from 3SAT to k-Clique.
 * The parameter k' = sat->num_clauses.
 * Runtime: O(m^2 * 7^2) = O(m^2) where m = num_clauses.
 *
 * @param sat         3SAT instance
 * @param[out] new_k   clique parameter (= sat->num_clauses)
 * @return            constructed graph, or NULL on failure
 */
graph_t *reduce_3sat_to_kclique_fpt(const sat3_instance_t *sat,
                                     int32_t *new_k);

/**
 * Verify the correctness of the 3SAT -> k-Clique reduction.
 * Checks:
 * 1. No edges between vertices of the same clause group (enforced by construction)
 * 2. Edges only between compatible assignments
 * 3. The construction graph has size O(m * 7)
 *
 * @param sat      original 3SAT instance
 * @param g        constructed graph
 * @param k        clique parameter
 * @return         true if reduction properties hold
 */
bool verify_3sat_to_kclique_reduction(const sat3_instance_t *sat,
                                       const graph_t *g, int32_t k);

/**
 * Recover a satisfying assignment of the 3SAT formula from
 * a k-clique found in the constructed graph.
 *
 * @param sat       original 3SAT instance
 * @param clique    k-clique in the constructed graph
 * @param[out] assignment  recovered Boolean assignment (caller allocates
 *                         sat->num_variables bools)
 * @return          true if assignment is valid and satisfies sat
 */
bool recover_sat_assignment(const sat3_instance_t *sat,
                             const clique_t *clique,
                             bool *assignment);

/**
 * Free a 3SAT instance.
 */
void sat3_instance_free(sat3_instance_t *sat);

/* ================================================================
 * FPT Reduction Framework
 * ================================================================ */

/**
 * Verify that a given reduction is a valid FPT reduction.
 *
 * Checks three properties:
 * 1. Correctness: (x,k) in A iff (x',k') in B
 * 2. Parameter boundedness: k' <= g(k) for some computable g
 * 3. Runtime: the reduction runs in FPT time f(k)*|x|^{O(1)}
 *
 * @param reduction  the FPT reduction to verify
 * @return           true if all FPT properties are satisfied
 */
bool verify_fpt_reduction(const fpt_reduction_t *reduction);

/**
 * Compute the parameter blowup function g(k) for the standard
 * k-Clique W[1]-hardness reduction from Weighted-3SAT.
 *
 * @param k  parameter of source problem
 * @return   parameter k' in the target k-Clique instance
 */
int32_t fpt_param_blowup_kclique(int32_t k);

/**
 * Check whether a reduction from problem A to k-Clique preserves
 * (or improves) the known lower bound.
 *
 * If A requires n^{alpha} time under ETH and the reduction
 * blows up the instance, the implied k-Clique bound may weaken.
 * This function computes the net effect.
 *
 * @param source_alpha   lower bound exponent for source problem
 * @param reduction_time_exponent  exponent c in O(n^c) reduction time
 * @param param_blowup   parameter growth in reduction
 * @param[out] clique_alpha  implied lower bound exponent for k-Clique
 * @return               true if a meaningful bound is obtained
 */
bool compute_transitive_lower_bound(double source_alpha,
                                     double reduction_time_exponent,
                                     int32_t param_blowup,
                                     double *clique_alpha);

#ifdef __cplusplus
}
#endif

#endif /* KCLIQUE_REDUCTION_H */
