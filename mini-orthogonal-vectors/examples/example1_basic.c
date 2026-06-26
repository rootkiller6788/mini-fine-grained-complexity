#include "ov.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(void) {
    srand((unsigned int)time(NULL));
    printf("=== Example 1: Basic Orthogonal Vectors ===\n\n");
    ov_instance_t *inst = ov_create_random(5, 6, 0.3, true);
    ov_print_summary(inst);
    ov_result_t r = ov_brute_force(inst);
    printf("Brute force: ");
    ov_report_result(&r);
    ov_destroy(inst);
    return 0;
}
