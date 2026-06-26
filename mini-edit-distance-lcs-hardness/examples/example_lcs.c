/* Example: LCS and Sequence Alignment */
#include <stdio.h>
#include <stdlib.h>
#include "../include/lcs.h"
#include "../include/alignment.h"

int main(void) {
    printf("=== LCS and Alignment Examples ===\n\n");

    printf("1. LCS computation:\n");
    lcs_result_t *r = lcs_compute("AGGTAB", "GXTXAYB");
    printf("   LCS('AGGTAB','GXTXAYB') = '%s' (len=%d)\n", r->subsequence, r->length);
    lcs_result_destroy(r);

    printf("\n2. Hirschberg linear-space LCS:\n");
    r = lcs_hirschberg("AGGTAB", "GXTXAYB");
    printf("   Hirschberg LCS = '%s' (len=%d)\n", r->subsequence, r->length);
    lcs_result_destroy(r);

    printf("\n3. Longest Common Substring:\n");
    printf("   LCSub('abcdef','zbcde') = %d\n", longest_common_substring("abcdef","zbcde"));

    printf("\n4. Shortest Common Supersequence:\n");
    printf("   SCS('ab','ac') length = %d\n", shortest_common_supersequence_length("ab","ac"));

    printf("\n5. LCS of 3 strings:\n");
    printf("   LCS3('abcde','ace','abce') = %d\n", lcs3_length("abcde","ace","abce"));

    printf("\n6. Needleman-Wunsch global alignment:\n");
    alignment_scoring_t scoring = {1, -1, -2, -1, false};
    alignment_result_t *aln = align_needleman_wunsch("AGCT", "AGGT", &scoring);
    printf("   NW score: %d\n", aln->total_score);
    alignment_result_print(aln);
    alignment_result_destroy(aln);

    printf("\n7. Smith-Waterman local alignment:\n");
    aln = align_smith_waterman("ACACACTA", "AGCACACA", &scoring);
    printf("   SW score: %d\n", aln->total_score);
    alignment_result_print(aln);
    alignment_result_destroy(aln);

    printf("\n8. DNA reverse complement:\n");
    char *rc = dna_reverse_complement("ATGCATGC");
    printf("   rev_comp('ATGCATGC') = '%s'\n", rc);
    free(rc);

    printf("\n9. GC content:\n");
    printf("   GC('ATGC') = %.1f%%\n", dna_gc_content("ATGC") * 100.0);

    return 0;
}
