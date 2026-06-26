#ifndef SUBQUADRATIC_H
#define SUBQUADRATIC_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <math.h>

/*
 * subquadratic.h -- Subquadratic Equivalence Class
 *
 * The subquadratic equivalence class contains problems whose best known
 * algorithms run in O(n^2) time and for which no O(n^{2-epsilon})
 * algorithm is known for any epsilon > 0.
 *
 * Canonical problem: Orthogonal Vectors (OV)
 * OV Conjecture: No O(n^{2-epsilon}) algorithm for OV when d = omega(log n).
 *
 * Theorem (Williams 2005): SETH => OV conjecture.
 * Therefore, SETH implies quadratic lower bounds for all problems in
 * the subquadratic equivalence class.
 *
 * Key equivalence results:
 *   OV ~ Edit Distance (Backurs & Indyk 2016)
 *   OV ~ LCS (Abboud, Backurs, Williams 2015)
 *   OV ~ Frechet Distance (Bringmann 2014)
 *   OV ~ DTW (Abboud, Backurs, Williams 2015)
 *
 * L1: OV, Edit Distance, LCS, Frechet, DTW definitions
 * L2: Subquadratic reduction concept
 * L3: Boolean vectors, dynamic programming tables, polygonal curves
 * L4: SETH => OV conjecture => Edit Distance / LCS lower bounds
 * L5: DP algorithms for Edit Distance, LCS, DTW
 * L6: OV, Edit Distance, LCS, Frechet, DTW
 */

/* ---- L1: Orthogonal Vectors Instance ---- */

typedef struct {
    int32_t     n;
    int32_t     d;
    uint32_t   *A;
    uint32_t   *B;
    int32_t     d_words;
} ov_instance_t;

/* ---- L1: Edit Distance (Levenshtein) Instance ---- */

typedef struct {
    int32_t     n;
    int32_t     m;
    char       *x;
    char       *y;
    int32_t    *dp;
    int32_t     distance;
} edit_distance_instance_t;

/* ---- L1: LCS Instance ---- */

typedef struct {
    int32_t     n;
    int32_t     m;
    char       *x;
    char       *y;
    int32_t    *dp;
    int32_t     lcs_length;
    char       *lcs_string;
} lcs_instance_t;

/* ---- L1: Frechet Distance Instance ---- */

typedef struct {
    int32_t     n;
    int32_t     m;
    double     *Px;
    double     *Py;
    double     *Qx;
    double     *Qy;
    double      distance;
} frechet_instance_t;

/* ---- L1: DTW (Dynamic Time Warping) Instance ---- */

typedef struct {
    int32_t     n;
    int32_t     m;
    double     *x;
    double     *y;
    double     *dp;
    double      dtw_distance;
} dtw_instance_t;

/* ---- L2/L5: Orthogonal Vectors Algorithms ---- */

ov_instance_t *ov_create(int32_t n, int32_t d);
void ov_free(ov_instance_t *inst);
void ov_set_bit(ov_instance_t *inst, int32_t vec_idx, int32_t bit_pos, bool value, bool in_A);
bool ov_get_bit(const ov_instance_t *inst, int32_t vec_idx, int32_t bit_pos, bool in_A);

/*
 * OV brute force: check every pair (a in A, b in B).
 * Time: O(n^2 * d / wordsize).
 */
bool ov_brute_force(const ov_instance_t *inst, int32_t *a_idx, int32_t *b_idx);

/*
 * Williams' OV algorithm for small dimension.
 * O(n^{2 - 1/O(log d/log n)}) when d = O(log n).
 * Reference: Williams (2005), ICALP.
 */
bool ov_williams_algorithm(const ov_instance_t *inst, int32_t *a_idx, int32_t *b_idx);

/*
 * Generate random OV instance with given density.
 */
void ov_random_instance(ov_instance_t *inst, double density, uint64_t seed);

/* ---- L2/L5: Edit Distance Algorithms ---- */

edit_distance_instance_t *edit_distance_create(int32_t n, int32_t m);
void edit_distance_free(edit_distance_instance_t *inst);

/*
 * Standard DP for edit distance (Levenshtein distance).
 * Operations: insert, delete, substitute (each cost 1).
 * Time: O(n*m), Space: O(n*m).
 *
 * Theorem (Backurs & Indyk 2016): Under SETH, no O(n^{2-epsilon})
 * algorithm for edit distance.
 */
int32_t edit_distance_dp(edit_distance_instance_t *inst);

/*
 * Edit distance with O(min(n,m)) space.
 * Two-row DP optimization.
 * Time: O(n*m), Space: O(min(n,m)).
 */
int32_t edit_distance_linear_space(const edit_distance_instance_t *inst);

/*
 * Simple Levenshtein distance on C strings.
 */
int32_t edit_distance_levenshtein(const char *x, int32_t n, const char *y, int32_t m);

/*
 * Edit distance with custom operation costs.
 */
int32_t edit_distance_weighted(const edit_distance_instance_t *inst,
                                int cost_ins, int cost_del, int cost_sub);

/* ---- L2/L5: LCS Algorithms ---- */

lcs_instance_t *lcs_create(int32_t n, int32_t m);
void lcs_free(lcs_instance_t *inst);

/*
 * Standard DP for LCS. Time: O(n*m), Space: O(n*m).
 *
 * Theorem (Abboud, Backurs, Williams 2015): Under SETH, no O(n^{2-epsilon})
 * algorithm for LCS, even for binary alphabet.
 */
int32_t lcs_dp(lcs_instance_t *inst);

/*
 * LCS with Hirschberg's algorithm: O(n*m) time, O(min(n,m)) space.
 * Reconstructs one actual LCS string.
 */
int32_t lcs_hirschberg(lcs_instance_t *inst);

/* ---- L2/L5: Frechet Distance Algorithms ---- */

frechet_instance_t *frechet_create(int32_t n, int32_t m);
void frechet_free(frechet_instance_t *inst);

/*
 * Continuous Frechet distance between two polygonal curves.
 * Uses the free-space diagram approach (Alt & Godau 1995).
 * Time: O(n*m).
 *
 * Theorem (Bringmann 2014): Under SETH, no O(n^{2-epsilon}) algorithm.
 */
double frechet_continuous(const frechet_instance_t *inst);

/* ---- L2/L5: DTW Algorithms ---- */

dtw_instance_t *dtw_create(int32_t n, int32_t m);
void dtw_free(dtw_instance_t *inst);

/*
 * Standard DP for DTW. Time: O(n*m).
 *
 * Theorem (Abboud, Backurs, Williams 2015): Under SETH, no O(n^{2-epsilon})
 * algorithm for DTW.
 */
double dtw_standard(dtw_instance_t *inst);

/* ---- L7: Applications ---- */

/*
 * DNA sequence alignment using edit distance.
 * Simulates Needleman-Wunsch global alignment.
 */
int32_t dna_sequence_alignment(const char *seq1, int32_t n1,
                                const char *seq2, int32_t n2,
                                int match_score, int mismatch_penalty,
                                int gap_penalty);

/*
 * Time series similarity search using DTW.
 */
int32_t time_series_search(const double *query, int32_t q_len,
                           const double **database, int32_t k, int32_t n,
                           double *best_distance);

/* ---- L6/L8: Subquadratic Completeness ---- */

/*
 * Verify subquadratic completeness for a problem P.
 */
bool verify_subquadratic_completeness(int32_t n, int32_t d,
    bool (*solve_P)(void *instance, int32_t n),
    void *(*reduce_OV_to_P)(const ov_instance_t *ov));

/*
 * Print the subquadratic equivalence class summary.
 */
void subquadratic_status_report(void);

#endif /* SUBQUADRATIC_H */