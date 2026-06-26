#include "ov.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(void) {
    srand((unsigned int)time(NULL));
    printf("=== OV Benchmark ===\n\n");
    int32_t ns[] = {50, 100, 200};
    int32_t d = 30;
    printf("Benchmarking OV with d=%d, density=0.3\n\n", d);
    for (int i = 0; i < 3; i++) {
        int32_t n = ns[i];
        printf("n=%d: ", n);
        ov_instance_t *inst = ov_create_random(n, d, 0.3, false);
        if (!inst) { printf("FAIL\n"); continue; }
        ov_result_t r = ov_benchmark(inst, ov_brute_force, 3);
        printf("brute=%.4fs, williams=", r.time_s);
        r = ov_benchmark(inst, ov_williams, 3);
        printf("%.4fs, ", r.time_s);
        r = ov_benchmark(inst, ov_four_russians, 3);
        printf("4russ=%.4fs\n", r.time_s);
        ov_destroy(inst);
    }
    return 0;
}
