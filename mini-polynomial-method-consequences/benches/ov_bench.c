#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "poly_method.h"

/* ============================================================================
 * ov_bench.c - OV Solver Benchmarks
 * ============================================================================ */

static double get_time_ms(void) {
    return (double)clock() * 1000.0 / (double)CLOCKS_PER_SEC;
}

int main(void) {
    printf("=== OV Solver Benchmarks ===\n\n");

    int32_t d = 16;

    for (int32_t n = 50; n <= 400; n *= 2) {
        ov_instance_t *A = ov_create(n, d);
        ov_instance_t *B = ov_create(n, d);

        uint64_t state = 12345;
        uint64_t mask = (1ULL << d) - 1;
        for (int32_t i = 0; i < n; i++) {
            state = state * 6364136223846793005ULL + 1442695040888963407ULL;
            A->vectors[i] = state & mask;
            state = state * 6364136223846793005ULL + 1442695040888963407ULL;
            B->vectors[i] = state & mask;
        }

        /* Packed method */
        double t1 = get_time_ms();
        ov_result_t r1 = ov_solve_packed(A, B);
        double t2 = get_time_ms();
        printf("n=%-4d packed:  %.3f ms  (found=%d)\n",
               n, t2 - t1, r1.orthogonal_pair_exists);

        /* Brute force */
        double t3 = get_time_ms();
        ov_result_t r2 = ov_solve_brute_force(A, B);
        double t4 = get_time_ms();
        printf("n=%-4d brute:   %.3f ms  (found=%d)\n",
               n, t4 - t3, r2.orthogonal_pair_exists);

        /* Count orthogonal pairs */
        double t5 = get_time_ms();
        int64_t cnt = ov_count_orthogonal_pairs(A, B);
        double t6 = get_time_ms();
        printf("n=%-4d count:   %.3f ms  (pairs=%lld)\n",
               n, t6 - t5, (long long)cnt);

        ov_destroy(A);
        ov_destroy(B);
    }

    printf("\n=== Benchmarks complete ===\n");
    return 0;
}
