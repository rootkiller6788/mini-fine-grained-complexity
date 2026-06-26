#include "ov.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(void) {
    srand((unsigned int)time(NULL));
    printf("=== Example 3: OV Applications ===\n\n");

    /* Database anti-join */
    ov_instance_t *db = ov_create_random(8, 6, 0.3, true);
    int32_t ra[20], rb[20];
    int64_t n = ov_database_anti_join(db, ra, rb, 20);
    printf("Database anti-join: %lld disjoint row pairs found\n", (long long)n);

    /* CRISPR off-target */
    const char *dna[] = {"ATCGATCGAT", "GCGCGCGCGC", "ATATATATAT"};
    ov_instance_t *crispr = ov_from_dna_sequences(dna, 3, 10);
    int32_t safe = ov_count_safe_crispr_pairs(crispr);
    printf("CRISPR safe pairs: %d\n", safe);
    ov_destroy(crispr);

    /* Code dissimilarity */
    int32_t fa, fb;
    double dissim = ov_max_code_dissimilarity(db, &fa, &fb);
    printf("Max code dissimilarity: %.3f (fragments %d, %d)\n", dissim, fa, fb);

    ov_destroy(db);
    return 0;
}
