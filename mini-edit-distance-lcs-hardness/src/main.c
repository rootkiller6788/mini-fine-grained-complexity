/* ============================================================================
 * main.c -- Entry point for mini-edit-distance-lcs-hardness
 *
 * Demonstrates edit distance, LCS, alignment, Frechet distance,
 * and SETH conditional lower bounds.
 * ============================================================================ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "edit_distance.h"
#include "lcs.h"
#include "seth_hardness.h"
#include "alignment.h"
#include "frechet_distance.h"
#include "string_utils.h"

int main(void) {
    printf("=== mini-edit-distance-lcs-hardness ===\n\n");

    /* Edit Distance Demo */
    printf("--- Edit Distance ---\n");
    const char *s1 = "kitten";
    const char *s2 = "sitting";
    int32_t ed = edit_distance_full(s1, s2);
    printf("ED('%s','%s') = %d\n", s1, s2, ed);

    ed = edit_distance_linear_space(s1, s2);
    printf("ED_linear('%s','%s') = %d\n", s1, s2, ed);

    ed = damerau_levenshtein("abc", "acb");
    printf("DL('abc','acb') = %d\n", ed);

    double js = jaro_winkler("martha", "marhta", 0.1);
    printf("Jaro-Winkler('martha','marhta') = %.3f\n", js);

    /* LCS Demo */
    printf("\n--- LCS ---\n");
    int32_t lcs_len = lcs_length("AGGTAB", "GXTXAYB");
    printf("LCS_length('AGGTAB','GXTXAYB') = %d\n", lcs_len);

    lcs_result_t *lcs_r = lcs_compute("AGGTAB", "GXTXAYB");
    if (lcs_r) {
        printf("LCS = '%s' (len=%d)\n", lcs_r->subsequence, lcs_r->length);
        lcs_result_destroy(lcs_r);
    }

    /* Hirschberg Demo */
    lcs_result_t *hirsch = lcs_hirschberg("AGGTAB", "GXTXAYB");
    if (hirsch) {
        printf("Hirschberg LCS = '%s'\n", hirsch->subsequence);
        lcs_result_destroy(hirsch);
    }

    /* Alignment Demo */
    printf("\n--- Sequence Alignment ---\n");
    alignment_scoring_t scoring = {1, -1, -2, -1, false};
    alignment_result_t *aln = align_needleman_wunsch("AGCT", "AGGT", &scoring);
    if (aln) {
        printf("NW alignment score: %d\n", aln->total_score);
        alignment_result_destroy(aln);
    }

    /* Frechet Distance Demo */
    printf("\n--- Frechet Distance ---\n");
    double px[] = {0, 1, 2, 3};
    double py[] = {0, 1, 0, 1};
    double qx[] = {0, 0.5, 1.5, 2.5};
    double qy[] = {0, 0.5, 0.5, 0};
    curve_t *cp = curve_from_arrays(px, py, 4);
    curve_t *cq = curve_from_arrays(qx, qy, 4);
    if (cp && cq) {
        double frech = discrete_frechet_distance(cp, cq);
        printf("Discrete Frechet = %.4f\n", frech);
        double dtw = dtw_distance(cp, cq);
        printf("DTW distance = %.4f\n", dtw);
        curve_destroy(cp); curve_destroy(cq);
    }

    /* SETH Hardness Report */
    printf("\n--- SETH Lower Bounds ---\n");
    print_comprehensive_lower_bound_report();

    /* String Utilities Demo */
    printf("\n--- String Utilities ---\n");
    char *rand_dna = dna_random(20, 12345);
    printf("Random DNA: %s\n", rand_dna);
    printf("GC content: %.1f%%\n", dna_gc_content(rand_dna) * 100.0);
    free(rand_dna);

    printf("\n=== All demonstrations complete ===\n");
    return 0;
}
