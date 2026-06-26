#include "ov.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(void) {
    srand((unsigned int)time(NULL));
    printf("=== Example 2: Algorithm Comparison ===\n\n");
    ov_instance_t *inst = ov_create_random(10, 10, 0.25, false);
    ov_print_summary(inst);
    ov_benchmark_all(inst);
    printf("\nDot product statistics:\n");
    ov_result_t r = ov_brute_force(inst);
    printf("Min dot: %d\n", ov_min_distance_to_orthogonal(inst));
    printf("Max dot: %d\n", ov_max_overlap(inst));
    ov_destroy(inst);
    return 0;
}
