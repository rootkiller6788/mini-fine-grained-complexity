/* ov_reduction.h -- Orthogonal Vectors Hypothesis and Reductions
 *
 * OV Hypothesis: For all epsilon > 0, Orthogonal Vectors on n vectors in
 * {0,1}^d with d = omega(log n) requires n^{2-epsilon} time in the
 * word-RAM model with O(log n) bit words.
 *
 * This is one of the central hypotheses in fine-grained complexity.
 * It is implied by SETH and implies hardness for many fundamental
 * problems in string algorithms, graph theory, and computational geometry.
 *
 * Key reductions:
 *   SETH → OV        (Williams, 2005)
 *   OV → Edit Distance   (Backurs-Indyk, 2015)
 *   OV → LCS             (Abboud-Backurs-Williams, 2015)
 *   OV → Diameter        (Roditty-Williams, 2013)
 *   OV → Subgraph Isomorphism  (bringmann-gronlund-larson, 2017)
 *
 * References:
 *   Williams (2005): "A new algorithm for optimal 2-constraint satisfaction"
 *   Backurs & Indyk (2015): "Edit Distance Cannot Be Computed in
 *        Strongly Subquadratic Time (unless SETH is false)"
 *   Abboud, Backurs, Williams (2015): "Tight Hardness Results for LCS
 *        and Other Sequence Similarity Measures"
 */

#ifndef OV_REDUCTION_H
#define OV_REDUCTION_H

#include "condlb.h"
#include <stdint.h>

/* ============================================================================
 * L1: Definitions — Orthogonal Vectors Problem
 * ============================================================================ */

/* A d-dimensional binary vector */
typedef struct {
    int   dim;       /* dimension d */
    int*  bits;      /* array of d bits (0 or 1) */
    int   id;        /* unique identifier */
} BinaryVector;

/* A set of n binary vectors in {0,1}^d */
typedef struct {
    BinaryVector* vectors;
    int           n_vectors;
    int           dim;
    int           capacity;
} VectorSet;

/* Result of checking orthogonality: are two vectors orthogonal?
 * Two vectors u,v in {0,1}^d are orthogonal iff their dot product is 0,
 * i.e., there is no coordinate i where both u_i = 1 and v_i = 1. */
typedef struct {
    int      i;          /* index of first vector */
    int      j;          /* index of second vector */
    int      orthogonal; /* 1 if orthogonal, 0 otherwise */
    uint64_t dot_product; /* sum_i u_i * v_i */
    int      found;      /* 1 if a solution was found */
} OrthogonalityResult;

/* ============================================================================
 * L5: Algorithms — OV Detection and Hardness Testing
 * ============================================================================ */

/* Brute-force OV: check all n choose 2 pairs. O(n^2 * d) time.
 * This is the baseline algorithm that the OV hypothesis claims is optimal. */
int ov_brute_force(const VectorSet* set, OrthogonalityResult* result);

/* Split-and-list OV for d = c*log n: partition dimensions, enumerate
 * all assignments. Time: O~(n^{2 - O(1/log c)}). Williams (2005).
 * This shows that OV becomes easier when d = O(log n). */
int ov_split_and_list(const VectorSet* set, double c_log_factor,
                      OrthogonalityResult* result);

/* Light-bulb OV: for sparse vectors (few 1s), use collision detection.
 * Time: O(n^{2-epsilon}) when density is low. */
int ov_light_bulb(const VectorSet* set, OrthogonalityResult* result);

/* FFT-based OV for moderate dimensions: use fast convolution.
 * Time: O~(n * 2^{d/2}). Advantageous when d < 2 log n. */
int ov_fft_method(const VectorSet* set, OrthogonalityResult* result);

/* ============================================================================
 * L5: OV Reductions — From OV to Other Problems
 * ============================================================================ */

/* OV → Edit Distance reduction (Backurs-Indyk, 2015).
 * Given OV instance (A, B, d), construct strings S1, S2 such that
 * EDIT(S1, S2) is small iff A and B have an orthogonal pair.
 * This shows: if EDIT is in O(n^{2-epsilon}), then OV is in O(n^{2-epsilon'}). */
void ov_to_edit_distance(const VectorSet* set,
                         char** out_s1, int* len1,
                         char** out_s2, int* len2);

/* OV → LCS reduction (Abboud-Backurs-Williams, 2015).
 * Constructs two strings whose LCS length reveals orthogonality. */
void ov_to_lcs(const VectorSet* set,
               char** out_s1, int* len1,
               char** out_s2, int* len2);

/* OV → Diameter reduction (Roditty-Williams, 2013).
 * Given OV instance, construct a graph whose diameter reveals orthogonality. */
void ov_to_diameter(const VectorSet* set,
                    int** out_adj_matrix, int* n_vertices);

/* OV → Subgraph Isomorphism (bringmann-gronlund-larson, 2017). */
void ov_to_subgraph_isomorphism(const VectorSet* set,
                                int** out_host_adj, int* host_n,
                                int** out_pat_adj,  int* pat_n);

/* ============================================================================
 * L2: Core Concepts — Vector Set Operations
 * ============================================================================ */

/* Create/destroy vector sets */
VectorSet* vs_create(int capacity, int dim);
void       vs_free(VectorSet* set);

/* Add a vector to the set. Returns its index. */
int vs_add_vector(VectorSet* set, const int* bits);

/* Generate random vector set of size n with specified dimension and density p */
VectorSet* vs_random(int n, int dim, double density);

/* Compute the dot product of two vectors in the set. */
uint64_t vs_dot_product(const VectorSet* set, int i, int j);

/* Check if any orthogonal pair exists (existence version of OV). */
int vs_has_orthogonal_pair(const VectorSet* set);

/* Find all orthogonal pairs. Returns number of pairs found. */
int vs_find_all_pairs(const VectorSet* set, OrthogonalityResult* results, int max_results);

/* ============================================================================
 * L7: Applications — SETH-based Lower Bounds
 * ============================================================================ */

/* Dynamic string matching: OV implies quadratic lower bound for
 * updates/query trade-offs in string matching. */
int ov_dynamic_string_matching_bound(const VectorSet* set,
                                      int* op_count, int* query_count);

/* Graph reachability oracle: OV implies hardness of constructing
 * small reachability oracles with fast query time. */
int ov_reachability_oracle_bound(const VectorSet* set,
                                  double* space_bound,
                                  double* query_time_bound);

/* Regular expression matching: OV implies quadratic lower bound
 * for regex matching (Backurs-Indyk, 2016). */
int ov_regex_matching_bound(const VectorSet* set,
                             double* lower_bound_exponent);

#endif /* OV_REDUCTION_H */
