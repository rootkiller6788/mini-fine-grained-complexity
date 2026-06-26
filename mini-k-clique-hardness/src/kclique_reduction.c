/**
 * kclique_reduction.c ? Reductions involving k-Clique
 *
 * Classical and parameterized reductions:
 * - k-Clique <-> k-Independent Set (complement)
 * - k-Clique -> k-Vertex Cover
 * - 3SAT -> k-Clique (FPT, W[1]-hardness proof)
 *
 * Reference: Karp, "Reducibility Among Combinatorial Problems" (1972)
 * Reference: Downey & Fellows, "Parameterized Complexity" (1999), Ch. 7
 * Reference: Flum & Grohe, "Parameterized Complexity Theory" (2006), Ch. 5
 *
 * Knowledge: L4 Fundamental Laws ? W[1]-completeness via reduction
 * Knowledge: L5 Algorithms ? Reduction constructions
 * Knowledge: L6 Canonical Problems ? k-IS, k-VC, 3SAT reductions
 */

#include "kclique_reduction.h"
#include "kclique_core.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ================================================================
 * k-Clique <-> k-Independent Set
 *
 * Theorem (folklore, implicit in Karp 1972):
 * G has a k-clique iff complement(G) has a k-independent set.
 *
 * This is parameter-preserving (k' = k) and runs in O(n^2) time.
 * Since k-Independent Set is W[1]-complete, so is k-Clique.
 * ================================================================ */

graph_t *kclique_to_kindependentset(const graph_t *g) {
    return graph_complement(g);
}

graph_t *kindependentset_to_kclique(const graph_t *g) {
    /* The complement of the complement is the original.
       So k-IS on G is equivalent to k-Clique on complement(G).
       This function is identical to kclique_to_kindependentset
       because both problems use the same reduction. */
    return graph_complement(g);
}

/* ================================================================
 * k-Clique -> k-Vertex Cover
 *
 * Theorem: G has a k-clique iff complement(G) has a vertex cover
 * of size n - k.
 *
 * Proof:
 * Let C be a k-clique in G. Then V\C is a vertex cover in
 * complement(G) because any edge in complement(G) has at least
 * one endpoint outside C (otherwise the edge would be in G[sic]).
 * Actually: C is a clique in G, so C is an independent set in
 * complement(G). Therefore V\C is a vertex cover in complement(G).
 * Size: |V| - |C| = n - k.
 * ================================================================ */

graph_t *kclique_to_vertexcover(const graph_t *g, int32_t k, int32_t *vc_k) {
    if (!g || !vc_k || k < 0 || k > g->n) return NULL;

    graph_t *comp = graph_complement(g);
    if (!comp) return NULL;

    *vc_k = g->n - k; /* vertex cover parameter */
    return comp;
}

/**
 * Reverse reduction: a vertex cover of size k' in G corresponds to
 * a clique of size n - k' in complement(G).
 *
 * @param g    original graph (for Vertex Cover)
 * @param vc_k vertex cover parameter
 * @param[out] clique_k  clique parameter in complement
 * @return     complement graph for clique search
 */
graph_t *vertexcover_to_kclique(const graph_t *g, int32_t vc_k, int32_t *clique_k) {
    if (!g || !clique_k || vc_k < 0 || vc_k > g->n) return NULL;

    graph_t *comp = graph_complement(g);
    if (!comp) return NULL;

    *clique_k = g->n - vc_k;
    return comp;
}

/* ================================================================
 * 3SAT -> k-Clique FPT Reduction (Canonical W[1]-hardness)
 *
 * This is THE reduction that proves k-Clique is W[1]-hard.
 *
 * Construction (Downey & Fellows 1995):
 *
 * Input: 3SAT formula phi with m clauses, each with 3 literals.
 *
 * Graph construction:
 * - For each clause C_j, create a "column" of 7 vertices, one per
 *   satisfying partial assignment to the 3 literals in C_j.
 *   (A partial assignment sets the 3 variables to make the clause true.)
 *
 * - Vertices: v_{j, a} where j indexes the clause (0..m-1) and
 *   a indexes the 7 satisfying assignments (0..6).
 *
 * - Edge (v_{j,a}, v_{j',a'}) exists iff:
 *   1. j != j' (different clauses)
 *   2. The partial assignments a and a' are compatible:
 *      they do not assign contradictory values to any variable.
 *
 * Claim: phi is satisfiable iff the constructed graph has an m-clique.
 *
 * Proof:
 * => If phi is satisfiable, choose the assignment for each clause
 *    that agrees with the global satisfying assignment. These m vertices
 *    are from different clauses and pairwise compatible => m-clique.
 *
 * <= An m-clique has exactly one vertex per clause (since vertices in
 *    the same clause are not adjacent). The m partial assignments are
 *    pairwise compatible, so their union gives a satisfying assignment.
 *
 * Parameter: k = m (number of clauses), so this is a parameterized
 * reduction with k' = m.
 *
 * Runtime: O(m^2 * 49) = O(m^2), which is polynomial.
 * ================================================================ */

/* Literal encoding: positive var i => i, negative var i => -(i+1) */
static int32_t __attribute__((unused)) encode_literal(int32_t var, bool positive) {
    return positive ? var : -(var + 1);
}

static int32_t decode_var(int32_t lit) {
    return (lit >= 0) ? lit : (-lit - 1);
}

static bool decode_sign(int32_t lit) {
    return lit >= 0;
}

/* Generate all 7 satisfying assignments for a 3-literal clause */
static void clause_satisfying_assignments(const int32_t lits[3],
                                           int32_t assignments[7][3]) {
    int idx = 0;
    for (int mask = 1; mask < 8; mask++) {
        /* mask bits: which literals are set to true */
        for (int b = 0; b < 3; b++) {
            assignments[idx][b] = (mask & (1 << b)) ? lits[b] : 0;
        }
        idx++;
    }
}

static bool assignments_compatible(const int32_t a1[3], const int32_t a2[3]) {
    for (int i = 0; i < 3; i++) {
        int32_t lit1 = a1[i];
        if (lit1 == 0) continue; /* not set in this assignment */
        for (int j = 0; j < 3; j++) {
            int32_t lit2 = a2[j];
            if (lit2 == 0) continue;
            int32_t v1 = decode_var(lit1), v2 = decode_var(lit2);
            if (v1 == v2) {
                /* Same variable: must have same sign */
                if (decode_sign(lit1) != decode_sign(lit2)) return false;
            }
        }
    }
    return true;
}

graph_t *reduce_3sat_to_kclique_fpt(const sat3_instance_t *sat,
                                     int32_t *new_k) {
    if (!sat || !new_k || sat->num_clauses <= 0 || sat->num_variables <= 0) {
        if (new_k) *new_k = 0;
        return NULL;
    }

    int32_t m = sat->num_clauses;
    *new_k = m;

    /* Total vertices: 7 * m (7 per clause) */
    int32_t n = 7 * m;
    graph_t *g = graph_create(n);
    if (!g) return NULL;

    /* Generate assignment tables for all clauses */
    int32_t (*assignments)[7][3] = (int32_t (*)[7][3])
        malloc((size_t)m * 7 * 3 * sizeof(int32_t));
    if (!assignments) { graph_destroy(g); return NULL; }

    for (int32_t c = 0; c < m; c++) {
        int32_t lits[3];
        for (int32_t i = 0; i < 3; i++) {
            lits[i] = sat->clause_literals[3 * c + i];
        }
        clause_satisfying_assignments(lits, assignments[c]);
    }

    /* Add edges between compatible assignments from different clauses */
    for (int32_t c1 = 0; c1 < m; c1++) {
        for (int32_t a1 = 0; a1 < 7; a1++) {
            int32_t v1 = c1 * 7 + a1;
            for (int32_t c2 = c1 + 1; c2 < m; c2++) {
                for (int32_t a2 = 0; a2 < 7; a2++) {
                    int32_t v2 = c2 * 7 + a2;
                    if (assignments_compatible(assignments[c1][a1],
                                                assignments[c2][a2])) {
                        graph_add_edge(g, v1, v2);
                    }
                }
            }
        }
    }

    free(assignments);
    return g;
}

bool verify_3sat_to_kclique_reduction(const sat3_instance_t *sat,
                                       const graph_t *g, int32_t k) {
    if (!sat || !g) return false;
    if (k != sat->num_clauses) return false;
    if (g->n != 7 * sat->num_clauses) return false;

    /* Check: no edges within same clause group (vertices 7c..7c+6) */
    for (int32_t c = 0; c < sat->num_clauses; c++) {
        for (int32_t a1 = 0; a1 < 7; a1++) {
            for (int32_t a2 = a1 + 1; a2 < 7; a2++) {
                int32_t v1 = c * 7 + a1, v2 = c * 7 + a2;
                if (graph_has_edge(g, v1, v2)) return false;
            }
        }
    }

    return true;
}

bool recover_sat_assignment(const sat3_instance_t *sat,
                             const clique_t *clique,
                             bool *assignment) {
    if (!sat || !clique || !assignment || !clique->found) return false;
    if (clique->k != sat->num_clauses) return false;

    /* Initialize all variables to false */
    for (int32_t v = 0; v < sat->num_variables; v++) {
        assignment[v] = false;
    }

    /* Process each vertex in the clique */
    for (int32_t c = 0; c < clique->k; c++) {
        int32_t vertex = clique->vertices[c];
        int32_t clause_idx = vertex / 7;
        /* Get the literals of this clause */
        int32_t lits[3];
        for (int32_t i = 0; i < 3; i++) {
            lits[i] = sat->clause_literals[3 * clause_idx + i];
        }
        /* Generate the assignment that this vertex represents */
        int32_t assignments[7][3];
        clause_satisfying_assignments(lits, assignments);
        int32_t a = vertex % 7;

        for (int32_t i = 0; i < 3; i++) {
            int32_t lit = assignments[a][i];
            if (lit == 0) continue;
            int32_t var = decode_var(lit);
            assignment[var] = decode_sign(lit);
        }
    }

    /* Verify that the assignment satisfies the formula */
    for (int32_t c = 0; c < sat->num_clauses; c++) {
        bool clause_sat = false;
        for (int32_t i = 0; i < 3; i++) {
            int32_t lit = sat->clause_literals[3 * c + i];
            int32_t var = decode_var(lit);
            bool val = decode_sign(lit);
            if (assignment[var] == val) {
                clause_sat = true;
                break;
            }
        }
        if (!clause_sat) return false;
    }

    return true;
}

void sat3_instance_free(sat3_instance_t *sat) {
    if (!sat) return;
    free(sat->clause_literals);
    sat->clause_literals = NULL;
    sat->num_variables = 0;
    sat->num_clauses = 0;
}

/* ================================================================
 * FPT Reduction Verification
 * ================================================================ */

bool verify_fpt_reduction(const fpt_reduction_t *reduction) {
    if (!reduction) return false;
    /* An FPT reduction must satisfy:
       1. is_fpt flag set
       2. Non-zero runtime exponent (polynomial factor)
       3. Finite parameter blowup */
    return reduction->is_fpt &&
           reduction->runtime_exponent > 0.0 &&
           reduction->param_blowup >= 0;
}

int32_t fpt_param_blowup_kclique(int32_t k) {
    /* Standard 3SAT -> k-Clique reduction:
       parameter k' = num_clauses. For 3SAT parameterized by
       variable count k, the number of clauses can be up to O(k^3).
       This gives k' = O(k^3). */
    if (k <= 0) return 0;
    return k * k * k; /* cubic blowup */
}

bool compute_transitive_lower_bound(double source_alpha,
                                     double reduction_time_exponent,
                                     int32_t param_blowup,
                                     double *clique_alpha) {
    if (!clique_alpha) return false;

    /* If problem A requires n^{source_alpha} time, and reduction
       from A to B takes n^{reduction_time_exponent} time producing
       instance of size n' with parameter k' = param_blowup,
       then B requires n'^{source_alpha / param_blowup - ...} time.

       Simplified: the lower bound transfers but may be weakened
       by the parameter blowup. */

    if (source_alpha <= 0.0) { *clique_alpha = 0.0; return false; }

    double transferred = source_alpha / (double)(param_blowup > 0 ? param_blowup : 1);
    /* The reduction itself adds overhead */
    transferred -= reduction_time_exponent;
    if (transferred < 0.0) transferred = 0.0;

    *clique_alpha = transferred;
    return transferred > 0.0;
}
