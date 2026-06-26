/* string_hardness.h -- Conditional Lower Bounds for String Problems
 *
 * String problems form a major class in fine-grained complexity.
 * Under SETH/OVC, many classic string problems require quadratic time:
 *
 * Edit Distance: minimum number of insertions, deletions, and
 *   substitutions to transform string A into string B.
 *   CLB: n^{2-o(1)} under SETH (Backurs-Indyk, 2015).
 *
 * Longest Common Subsequence (LCS): longest sequence appearing
 *   as a subsequence in both strings.
 *   CLB: n^{2-o(1)} under SETH (Abboud-Backurs-Williams, 2015).
 *
 * Frechet Distance: measure of similarity between curves.
 *   CLB: n^{2-o(1)} under SETH (Bringmann, 2014).
 *
 * Dynamic Time Warping (DTW): sequence alignment allowing
 *   non-linear warping.
 *   CLB: n^{2-o(1)} under SETH (Abboud-Backurs-Williams, 2015;
 *                                 Bringmann-Kunnemann, 2015).
 *
 * Longest Palindrome Subsequence, Regular Expression Matching,
 * Pattern Matching with Wildcards, etc.
 *
 * References:
 *   Backurs & Indyk (2015): "Edit Distance Cannot Be Computed in
 *        Strongly Subquadratic Time (unless SETH is false)"
 *   Bringmann (2014): "Why Walking the Dog Takes Time: Frechet
 *        Distance Has No Strongly Subquadratic Algorithms Unless SETH Fails"
 *   Abboud, Backurs, Williams (2015): "Tight Hardness Results for
 *        LCS and Other Sequence Similarity Measures"
 */

#ifndef STRING_HARDNESS_H
#define STRING_HARDNESS_H

#include "condlb.h"
#include <stdint.h>

/* ============================================================================
 * L1: Definitions - String Distance Measures
 * ============================================================================ */

/* Edit operation types */
typedef enum {
    EDIT_NONE   = 0,
    EDIT_INSERT = 1,
    EDIT_DELETE = 2,
    EDIT_SUBSTITUTE = 3
} EditOperation;

/* Single edit operation in an alignment */
typedef struct {
    EditOperation op;
    int           pos1;  /* position in string 1 */
    int           pos2;  /* position in string 2 */
    char          old_char;
    char          new_char;
} EditStep;

/* Edit distance result with alignment trace */
typedef struct {
    int        distance;
    int        n_steps;
    EditStep*  steps;
} EditResult;

/* LCS result */
typedef struct {
    int    length;
    char*  subsequence;
} LcsResult;

/* ============================================================================
 * L5: Algorithms - String Distance Computations
 * ============================================================================ */

/* Classic O(n*m) dynamic programming for edit distance.
 * Wagner-Fischer algorithm.
 * Time: O(n*m), Space: O(min(n,m)).
 * Reference: Wagner & Fischer (1974) */
int edit_distance_dp(const char* s1, int len1,
                     const char* s2, int len2);

/* Needleman-Wunsch global alignment (edit distance with alignment trace).
 * Time: O(n*m), Space: O(n*m).
 * Reference: Needleman & Wunsch (1970) */
EditResult edit_distance_align(const char* s1, int len1,
                               const char* s2, int len2);

/* Four-Russians speedup for edit distance: divide strings into blocks,
 * precompute all possible block edit distances.
 * Time: O(n^2 / log^2 n) when alphabet is small.
 * Reference: Masek & Paterson (1980) */
int edit_distance_four_russians(const char* s1, int len1,
                                const char* s2, int len2,
                                int block_size);

/* O(n + d^2) edit distance using Ukkonen's algorithm.
 * Where d is the actual edit distance.
 * Efficient when strings are similar.
 * Time: O(n + d^2), Space: O(d).
 * Reference: Ukkonen (1985) */
int edit_distance_ukkonen(const char* s1, int len1,
                          const char* s2, int len2);

/* Classic O(n*m) LCS dynamic programming.
 * Time: O(n*m), Space: O(min(n,m)). */
int lcs_length(const char* s1, int len1,
               const char* s2, int len2);

/* LCS with full subsequence extraction.
 * Time: O(n*m), Space: O(n*m). */
LcsResult lcs_extract(const char* s1, int len1,
                      const char* s2, int len2);

/* Hunt-Szymanski LCS for strings with few matches.
 * Time: O((R+n) log n) where R = number of matching pairs.
 * Reference: Hunt & Szymanski (1977) */
int lcs_hunt_szymanski(const char* s1, int len1,
                       const char* s2, int len2);

/* Frechet distance between two polygonal curves.
 * Discrete Frechet distance via DP.
 * Time: O(n*m).
 * Reference: Eiter & Mannila (1994) */
double frechet_distance(const double* px1, const double* py1, int n1,
                        const double* px2, const double* py2, int n2);

/* Dynamic Time Warping (DTW) distance.
 * Time: O(n*m).
 * Reference: Sakoe & Chiba (1978) */
double dtw_distance(const double* seq1, int n1,
                    const double* seq2, int n2);

/* ============================================================================
 * L4: Fundamental Theorems - Quadratic Lower Bounds
 * ============================================================================ */

/* Theorem (Backurs-Indyk, 2015):
 * If Edit Distance can be computed in O(n^{2-epsilon}) time,
 * then SETH is false (i.e., k-SAT can be solved in
 * O(2^{(1-delta)n}) for some delta > 0).
 *
 * This function computes the eps-delta tradeoff:
 * given epsilon (the subquadratic speedup for edit distance),
 * compute the implied delta (the SETH violation). */
double edit_distance_seth_tradeoff(double epsilon_speedup);

/* Theorem (Bringmann-Kunnemann, 2015):
 * Even strongly subquadratic approximation of Frechet distance
 * (within factor < 3) would refute SETH. */
int frechet_approx_bound(double approx_factor);

/* Theorem (Abboud-Backurs-Williams, 2015):
 * LCS requires n^{2-o(1)} time under SETH, even for
 * binary alphabet. */
double lcs_seth_tradeoff(double epsilon_speedup);

/* ============================================================================
 * L7: Applications
 * ============================================================================ */

/* Compute edit distance conditional lower bound for given input sizes.
 * Returns the exponent: edit_distance >= n^{return_value} under SETH. */
double string_clb_edit_distance(int n, int alphabet_size);

/* Check if a string algorithm's runtime is consistent with SETH.
 * Returns 1 if algorithm could refute SETH. */
int string_check_seth_consistency(const char* algo_name,
                                  double algo_exponent,
                                  int n);

/* Create a SETH-to-Edit-Distance hard instance:
 * Given a SAT instance, construct strings whose edit distance
 * would reveal satisfiability. */
void string_seth_to_edit_hard_instance(int n_vars, int n_clauses,
                                       char** out_s1, int* len1,
                                       char** out_s2, int* len2);

#endif /* STRING_HARDNESS_H */
