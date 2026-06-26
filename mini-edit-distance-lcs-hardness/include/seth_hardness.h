#ifndef SETH_HARDNESS_H
#define SETH_HARDNESS_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <math.h>

/* ============================================================================
 * mini-edit-distance-lcs-hardness: SETH Hardness Framework
 *
 * The Strong Exponential Time Hypothesis (SETH) framework applied to
 * edit distance and LCS. This header provides the conditional lower
 * bound structures and OV (Orthogonal Vectors) reduction machinery.
 *
 * SETH: ∀ε>0, ∃k s.t. k-SAT ∉ DTIME(2^{(1-ε)n})
 *
 * Key results:
 *   - Backurs-Indyk (STOC 2015): Edit Distance ∉ O(n^{2-ε}) unless SETH false
 *   - ABW (FOCS 2015): LCS ∉ O(n^{2-ε}) unless SETH false
 *   - Bringmann-Kunnemann (FOCS 2015): Frechet ∉ O(n^{2-ε}) unless SETH false
 *
 * L1: SETH, OV Conjecture definitions
 * L2: Fine-grained reduction framework
 * L3: Reduction graph structures
 * L4: Backurs-Indyk Theorem, ABW Theorem
 * L8: Advanced conditional lower bounds
 * ============================================================================ */

/* ---- L1: Core Conjecture Definitions ---- */

/* SETH parameters */
typedef struct {
    double   epsilon;      /* ε in SETH statement */
    int32_t  k;            /* k-SAT clause width */
    double   s_k;          /* s_k = inf{c: k-SAT in O(2^{cn})} */
    double   s_infinity;   /* s_∞ = lim_{k→∞} s_k = 1 */
} seth_statement_t;

/* Orthogonal Vectors (OV) Conjecture parameters */
typedef struct {
    int32_t  n;            /* Number of vectors */
    int32_t  d;            /* Dimension */
    double   target_time;  /* Claimed O(n^{2-ε}) time */
    double   epsilon;      /* ε > 0 */
} ov_conjecture_t;

/* ---- L2: OV Problem Data Structures ---- */

/* A binary vector */
typedef struct {
    int32_t  dim;          /* Dimension d */
    uint64_t *bits;        /* Packed bits: ceil(d/64) words */
    int32_t  num_words;
} binary_vector_t;

/* A set of binary vectors (one OV instance) */
typedef struct {
    binary_vector_t *set_a;
    int32_t          num_a;
    binary_vector_t *set_b;
    int32_t          num_b;
    int32_t          dim;
} ov_instance_t;

/* Orthogonal Vectors query result */
typedef struct {
    bool      exists_orthogonal;
    int32_t   index_a;
    int32_t   index_b;
    uint64_t  operations;   /* Number of comparisons performed */
} ov_result_t;

/* ---- OV Operations ---- */

/* Create and destroy binary vectors */
binary_vector_t *bv_create(int32_t dim);
void bv_destroy(binary_vector_t *v);
void bv_set_bit(binary_vector_t *v, int32_t pos, bool value);
bool bv_get_bit(const binary_vector_t *v, int32_t pos);
bool bv_are_orthogonal(const binary_vector_t *a, const binary_vector_t *b);
int32_t bv_dot_product(const binary_vector_t *a, const binary_vector_t *b);

/* OV instance operations */
ov_instance_t *ov_instance_create(int32_t num_a, int32_t num_b, int32_t dim);
void ov_instance_destroy(ov_instance_t *inst);
void ov_instance_random_fill(ov_instance_t *inst, uint64_t seed);

/* Brute-force OV: O(n^2 * d) time */
ov_result_t ov_brute_force(const ov_instance_t *inst);

/* OV with 4-list splitting: reduces OV(n,d) to 4 OV(n/2, d) subproblems */
ov_result_t ov_four_list_split(const ov_instance_t *inst);

/* ---- L3: Reduction Framework ---- */

/* Types of fine-grained reductions */
typedef enum {
    REDUCTION_OV_TO_EDIT = 0,
    REDUCTION_OV_TO_LCS  = 1,
    REDUCTION_SAT_TO_OV  = 2,
    REDUCTION_KCLIQUE_TO_LCS = 3,
    REDUCTION_3SUM_TO_EDIT   = 4
} reduction_type_t;

/* A fine-grained reduction with timing guarantees */
typedef struct {
    reduction_type_t type;
    int32_t          input_size;
    int32_t          output_size;
    double           input_time;    /* Time to solve input problem */
    double           output_time;   /* Time to solve output problem */
    double           reduction_time; /* Time to perform reduction */
    double           epsilon;        /* ε for n^{2-ε} claim */
    bool             is_tight;       /* Is the reduction tight? */
} reduction_info_t;

/* ---- SETH Conditional Lower Bound Verification ---- */

/* Check if an algorithm with given time complexity would refute SETH
 * for the edit distance problem.
 * Input: time_exponent c where algorithm runs in O(n^c)
 * Returns true if SETH would be refuted (i.e., c < 2) */
bool seth_edit_distance_check(double time_exponent);

/* Check if an algorithm would refute SETH for LCS */
bool seth_lcs_check(double time_exponent);

/* Compute the reduction overhead factor */
double reduction_overhead(int32_t input_n, int32_t output_n,
                          int32_t from_dim, int32_t to_dim);

/* ---- L4: Theorem Verification Functions ---- */

/* Backurs-Indyk (2015): Verify the reduction chain
 * SAT(n) → OV(N, d) → Edit Distance(N')
 *
 * Given N = 2^{n/2}, d = O(log N), verify that:
 *   Edit Distance in O(N^{2-ε}) ⇒ SAT in O(2^{(1-δ)n})
 *   for some δ > 0, contradicting SETH. */
bool backurs_indyk_verify(int32_t n, double epsilon);

/* ABW (2015): Verify the reduction chain
 * SAT(n) → OV(N, d) → LCS(N')
 *
 * Given N = 2^{n/2}, d = O(log N), verify that:
 *   LCS in O(N^{2-ε}) ⇒ SAT in O(2^{(1-δ)n})
 *   for some δ > 0, contradicting SETH. */
bool abw_verify(int32_t n, double epsilon);

/* Bringmann-Kunnemann (2015): Verify Frechet distance reduction
 *
 *   Frechet in O(N^{2-ε}) ⇒ SAT in O(2^{(1-δ)n}) */
bool bringmann_kunnemann_verify(int32_t n, double epsilon);

/* ---- L2: Conditional Lower Bound Reports ---- */

/* Print the SETH conditional lower bound status for edit distance */
void print_edit_seth_status(void);

/* Print the SETH conditional lower bound status for LCS */
void print_lcs_seth_status(void);

/* Generate a comprehensive lower bound report */
void print_comprehensive_lower_bound_report(void);

/* ---- L8: Beyond Quadratic Lower Bounds ---- */

/* Check if edit distance can be solved in truly sub-quadratic time
 * under the OV conjecture (without going through SETH). */
bool ov_to_edit_implies_subquadratic_edit(double time_exponent);

/* The k-OV problem: find k orthogonal vectors.
 * For k=2, this is the standard OV problem.
 * For k>2, this is the k-Orthogonal Vectors problem. */
bool k_ov_to_edit_implies(double time_exponent, int32_t k);

/* OVH (Orthogonal Vectors Hypothesis):
 * OV problem requires n^{2-o(1)} time for d = ω(log n). */
bool ovh_check(double time_exponent, int32_t d, int32_t n);

/* All-Pairs Shortest Paths (APSP) to Edit Distance reduction check */
bool apsp_to_edit_check(double time_exponent);

/* 3SUM to edit distance reduction check */
bool three_sum_to_edit_check(double time_exponent);

/* ---- L9: Research Frontier Utilities ---- */

/* Compute the "SETH constant" s_k for a given k.
 * s_k = inf{c : k-SAT ∈ DTIME(2^{c n})}
 * Currently, s_3 ≤ 0.3863 (PPSZ 2005).
 * SETH says: lim_{k→∞} s_k = 1. */
double seth_constant_s_k(int32_t k);

/* Check if an algorithm running in time O(n^c) would imply
 * breakthrough results for multiple problems */
void breakthrough_implications(double c);

#endif /* SETH_HARDNESS_H */
