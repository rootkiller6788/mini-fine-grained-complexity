/* Example: Subquadratic Equivalence -- OV, Edit Distance, LCS, DTW
 *
 * Demonstrates the subquadratic equivalence class by running:
 *   1. Orthogonal Vectors brute-force and Williams algorithm
 *   2. Edit Distance computation between two strings
 *   3. Longest Common Subsequence
 *   4. Dynamic Time Warping distance
 *
 * References:
 *   Backurs & Indyk (2016), STOC
 *   Abboud, Backurs, Williams (2015), FOCS
 *   Bringmann (2014), FOCS
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "subquadratic.h"

int main(void) {
    printf("=== Example: Subquadratic Equivalence Class ===\n\n");

    /* --- Orthogonal Vectors --- */
    printf("--- Orthogonal Vectors (OV) ---\n");
    int32_t n = 5, d = 4;
    ov_instance_t *ov = ov_create(n, d);

    /* Set up vectors: A has ones, B has complementary ones so no orthogonal pair */
    ov_set_bit(ov, 0, 0, true, true);
    ov_set_bit(ov, 0, 0, true, false);  /* Not orthogonal: both have bit 0 */

    ov_set_bit(ov, 1, 0, true, true);
    ov_set_bit(ov, 1, 1, true, false);  /* Orthogonal: bits don't overlap */

    int32_t ai, bi;
    bool found = ov_brute_force(ov, &ai, &bi);
    printf("Brute force OV: %s", found ? "found" : "not found");
    if (found) printf(" -- pair (A[%d], B[%d]) is orthogonal", ai, bi);
    printf("\n");

    ov_free(ov);

    /* --- Edit Distance --- */
    printf("\n--- Edit Distance ---\n");
    const char *s1 = "kitten";
    const char *s2 = "sitting";
    int32_t ed = edit_distance_levenshtein(s1, (int32_t)strlen(s1),
                                             s2, (int32_t)strlen(s2));
    printf("Edit Distance('%s', '%s') = %d\n", s1, s2, ed);
    printf("(Expected: 3 -- substitute k->s, e->i, insert g)\n");

    /* --- LCS --- */
    printf("\n--- Longest Common Subsequence ---\n");
    lcs_instance_t *lcs = lcs_create(7, 7);
    strncpy(lcs->x, "AGGTABX", 7);
    strncpy(lcs->y, "GXTXAYB", 7);
    int32_t lcs_len = lcs_hirschberg(lcs);
    printf("LCS('%s', '%s') = %d\n", lcs->x, lcs->y, lcs_len);
    if (lcs->lcs_string) printf("LCS string: %s\n", lcs->lcs_string);
    printf("(Expected: 4 -- GTAB)\n");
    lcs_free(lcs);

    /* --- DTW --- */
    printf("\n--- Dynamic Time Warping ---\n");
    dtw_instance_t *dtw = dtw_create(5, 5);
    double x[] = {1.0, 3.0, 4.0, 9.0, 2.0};
    double y[] = {2.0, 4.0, 5.0, 8.0, 1.0};
    for (int32_t i = 0; i < 5; i++) { dtw->x[i] = x[i]; dtw->y[i] = y[i]; }
    double dtw_dist = dtw_standard(dtw);
    printf("DTW distance: %.3f\n", dtw_dist);
    dtw_free(dtw);

    /* --- DNA Alignment --- */
    printf("\n--- DNA Sequence Alignment ---\n");
    int32_t score = dna_sequence_alignment("ACGTACGT", 8, "ACTACTGT", 8, 2, -1, -1);
    printf("DNA alignment score (ACGTACGT vs ACTACTGT): %d\n", score);

    /* --- Status Report --- */
    printf("\n--- Subquadratic Status ---\n");
    subquadratic_status_report();

    printf("\n=== Example Complete ===\n");
    return 0;
}