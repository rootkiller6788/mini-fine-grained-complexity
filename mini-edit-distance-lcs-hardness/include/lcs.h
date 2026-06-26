#ifndef LCS_H
#define LCS_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <math.h>

/* ============================================================================
 * mini-edit-distance-lcs-hardness: Longest Common Subsequence (LCS)
 *
 * Core definitions and data structures for LCS computation and
 * conditional lower bounds.
 *
 * Reference: Hirschberg (1975) "A linear space algorithm for computing
 *            maximal common subsequences" CACM 18(6):341-343
 *            Abboud, Backurs, Williams (FOCS 2015)
 *            "Tight Hardness Results for LCS and Other Sequence Similarity
 *            Measures"
 *
 * L1: Core definitions: LCS, subsequence, common subsequence
 * L2: Quadratic-time lower bound under SETH
 * L3: DP matrix for LCS, alignment representation
 * ============================================================================ */

/* ---- L1: Core Type Definitions ---- */

/* An aligned pair in the LCS */
typedef struct {
    int32_t pos_a;      /* Position in string A */
    int32_t pos_b;      /* Position in string B */
    char    character;  /* The matching character */
} lcs_match_t;

/* Result of LCS computation */
typedef struct {
    int32_t     length;           /* Length of the LCS */
    lcs_match_t *matches;         /* Array of matching positions */
    int32_t     num_matches;
    char       *subsequence;      /* The actual LCS string (null-terminated) */
} lcs_result_t;

/* ---- L2: DP-based LCS Computation ---- */

/* Standard DP: O(n*m) time, O(n*m) space.
 * Returns the length of the LCS. */
int32_t lcs_length(const char *a, const char *b);

/* Compute full LCS result including the actual subsequence.
 * O(n*m) time, O(n*m) space. */
lcs_result_t *lcs_compute(const char *a, const char *b);

/* Hirschberg's linear-space LCS algorithm.
 * O(n*m) time, O(min(n,m)) space.
 * Returns both length and the actual subsequence. */
lcs_result_t *lcs_hirschberg(const char *a, const char *b);

/* LCS length only, O(min(n,m)) space using 2-row DP.
 * Optimal for when only the length is needed. */
int32_t lcs_length_linear_space(const char *a, const char *b);

/* Free an LCS result */
void lcs_result_destroy(lcs_result_t *r);

/* Print the LCS result in human-readable format */
void lcs_result_print(const lcs_result_t *r);

/* ---- L3: LCS DP Table ---- */

typedef struct {
    int32_t *data;
    int32_t  rows;
    int32_t  cols;
} lcs_dp_table_t;

typedef struct {
    int32_t value;
    uint8_t direction;  /* 0=diagonal, 1=up, 2=left */
} lcs_dp_cell_t;

typedef struct {
    lcs_dp_cell_t *data;
    int32_t        rows;
    int32_t        cols;
} lcs_dp_traceable_t;

lcs_dp_table_t *lcs_dp_create(int32_t len_a, int32_t len_b);
void lcs_dp_destroy(lcs_dp_table_t *dp);
void lcs_dp_fill(lcs_dp_table_t *dp, const char *a, const char *b);

lcs_dp_traceable_t *lcs_dp_traceable_create(int32_t len_a, int32_t len_b);
void lcs_dp_traceable_destroy(lcs_dp_traceable_t *dp);
void lcs_dp_traceable_fill(lcs_dp_traceable_t *dp, const char *a, const char *b);

/* ---- Variants of LCS ---- */

/* Longest Common Substring (contiguous):
 * O(n*m) DP, returns length */
int32_t longest_common_substring(const char *a, const char *b);

/* Longest Common Substring with result */
typedef struct {
    int32_t  length;
    int32_t  start_a;
    int32_t  start_b;
    char    *substring;
} lcsubstr_result_t;

lcsubstr_result_t *lcsubstr_compute(const char *a, const char *b);
void lcsubstr_result_destroy(lcsubstr_result_t *r);

/* Shortest Common Supersequence (SCS) */
int32_t shortest_common_supersequence_length(const char *a, const char *b);
char *shortest_common_supersequence(const char *a, const char *b);

/* LCS of 3 strings: O(n*m*p) time */
int32_t lcs3_length(const char *a, const char *b, const char *c);

/* LCS of k strings via DP on position vectors */
int32_t lcs_k_length(char **strings, int32_t k);

/* ---- Sequence Similarity Measures ---- */

/* Edit distance via LCS (for unit-cost operations):
 * edit_dist = |a| + |b| - 2*|LCS(a,b)| */
int32_t edit_via_lcs(const char *a, const char *b);

/* Normalized LCS similarity: |LCS| / max(|a|, |b|) */
double lcs_similarity(const char *a, const char *b);

/* Longest increasing subsequence (LIS) in an integer array.
 * O(n log n) using patience sorting. */
int32_t longest_increasing_subsequence(const int32_t *arr, int32_t n);

/* Longest common increasing subsequence (LCIS) of two arrays */
int32_t longest_common_increasing_subsequence(const int32_t *a, int32_t n,
                                              const int32_t *b, int32_t m);

/* ---- L5: Faster LCS Algorithms ---- */

/* Hunt-Szymanski algorithm for LCS:
 * Efficient when alphabet size is large and matches are few.
 * O((R+n) log n) where R = number of matching pairs. */
lcs_result_t *lcs_hunt_szymanski(const char *a, const char *b);

/* Four Russians speedup applied to LCS:
 * Conceptual O(n^2 / log n) algorithm using block precomputation. */
int32_t lcs_four_russians(const char *a, const char *b);

/* ---- L8: LCS Hardness Framework ---- */

/* Verify if a proposed algorithm would refute SETH for edit/LCS */
bool lcs_would_refute_seth(double exponent, int32_t n);

/* Compute the alignment graph for OV-to-LCS reduction */
int32_t lcs_ov_reduction_check(int32_t n, int32_t d);

#endif /* LCS_H */
