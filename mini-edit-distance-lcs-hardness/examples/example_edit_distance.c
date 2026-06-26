/* Example: Edit Distance Computation */
#include <stdio.h>
#include <stdlib.h>
#include "../include/edit_distance.h"
#include "../include/string_utils.h"

int main(void) {
    printf("=== Edit Distance Examples ===\n\n");
    printf("1. Levenshtein distance:\n");
    printf("   ED('kitten','sitting') = %d\n", edit_distance_full("kitten","sitting"));
    printf("   ED('Saturday','Sunday') = %d\n", edit_distance_full("Saturday","Sunday"));
    printf("   ED('','abc') = %d\n", edit_distance_full("","abc"));

    printf("\n2. Hamming distance:\n");
    printf("   Hamming('abc','abd') = %d\n", hamming_distance("abc","abd"));

    printf("\n3. Damerau-Levenshtein:\n");
    printf("   DL('abc','acb') = %d (transposition)\n", damerau_levenshtein("abc","acb"));

    printf("\n4. Jaro-Winkler (record linkage):\n");
    printf("   JW('martha','marhta') = %.3f\n", jaro_winkler("martha","marhta",0.1));
    printf("   JW('Dwayne','Duane') = %.3f\n", jaro_winkler("Dwayne","Duane",0.1));

    printf("\n5. Edit script reconstruction:\n");
    edit_script_t *script = edit_distance_with_script("abc", "axbc");
    ed_script_print(script);
    ed_script_destroy(script);

    printf("\n6. Bounded edit distance (k=1):\n");
    printf("   ED_bounded('hello','hallo',1) = %d\n", edit_distance_bounded("hello","hallo",1));

    printf("\n7. Normalized edit distance:\n");
    printf("   NED('abc','abcd') = %.2f\n", edit_distance_normalized("abc","abcd"));

    printf("\n8. Random strings + edit distance:\n");
    char *a = str_random(20, "ACGT", 42);
    char *b = str_mutate(a, 0.2, 99);
    printf("   Original: %s\n", a);
    printf("   Mutated:  %s\n", b);
    printf("   ED = %d\n", edit_distance_full(a, b));
    free(a); free(b);

    return 0;
}
