/* Micro-benchmarks for equivalence class operations.
 * Measures runtime of key algorithms to compare with theoretical predictions. */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "subcubic.h"
#include "subquadratic.h"
#include "threesum.h"

static double now_sec(void) {
    return (double)clock() / CLOCKS_PER_SEC;
}

int main(void) {
    printf("=== Benchmarks: Equivalence Classes ===\n\n");

    /* APSP benchmark */
    int32_t n = 100;
    apsp_instance_t *apsp = apsp_create(n);
    srand(42);
    for (int32_t i = 0; i < n; i++)
        for (int32_t j = 0; j < n; j++)
            if (i != j && (rand() % 100) < 30)
                apsp_set_edge(apsp, i, j, (double)(rand() % 100));

    double t0 = now_sec();
    apsp_floyd_warshall(apsp);
    double t1 = now_sec();
    printf("APSP (Floyd-Warshall, n=%d): %.4f sec\n", n, t1 - t0);
    apsp_free(apsp);

    /* Edit Distance benchmark */
    int32_t str_len = 500;
    char *s1 = (char *)malloc((size_t)(str_len + 1));
    char *s2 = (char *)malloc((size_t)(str_len + 1));
    for (int32_t i = 0; i < str_len; i++) {
        s1[i] = 'A' + (rand() % 4);
        s2[i] = 'A' + (rand() % 4);
    }
    s1[str_len] = s2[str_len] = '\0';

    t0 = now_sec();
    int32_t ed = edit_distance_levenshtein(s1, str_len, s2, str_len);
    t1 = now_sec();
    printf("Edit Distance (n=m=%d): %.4f sec, result=%d\n", str_len, t1 - t0, ed);
    free(s1); free(s2);

    /* 3SUM benchmark */
    int32_t n3 = 500;
    threesum_instance_t *ts = threesum_create(n3);
    for (int32_t i = 0; i < n3; i++)
        ts->values[i] = (int64_t)(rand() % (n3 * n3)) - (n3 * n3 / 2);
    t0 = now_sec();
    threesum_quadratic(ts);
    t1 = now_sec();
    printf("3SUM (n=%d): %.4f sec\n", n3, t1 - t0);
    threesum_free(ts);

    printf("\n=== Benchmarks Complete ===\n");
    return 0;
}