#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "poly_method.h"

/* ============================================================================
 * williams_bench.c - Williams SAT Algorithm Benchmarks
 * ============================================================================ */

static double get_time_ms(void) {
    return (double)clock() * 1000.0 / (double)CLOCKS_PER_SEC;
}

int main(void) {
    printf("=== Williams SAT Algorithm Benchmarks ===\n\n");

    int32_t k = 3;

    for (int32_t n = 4; n <= 18; n += 2) {
        int32_t m = n * 2;  /* clause-to-variable ratio = 2 */

        williams_cnf_t *f = williams_generate_ksat(n, m, k, (uint64_t)n * 100);

        /* Brute force */
        double t1 = get_time_ms();
        williams_result_t r1 = williams_brute_force(f);
        double t2 = get_time_ms();

        /* Williams polynomial method */
        double t3 = get_time_ms();
        williams_result_t r2 = williams_solve(f);
        double t4 = get_time_ms();

        printf("n=%-3d m=%-3d:  brute=%.3f ms  williams=%.3f ms  match=%d\n",
               n, m, t2 - t1, t4 - t3, r1.satisfiable == r2.satisfiable);

        williams_cnf_destroy(f);
        free(r1.witness);
        free(r2.witness);
    }

    printf("\n=== Benchmarks complete ===\n");
    return 0;
}
