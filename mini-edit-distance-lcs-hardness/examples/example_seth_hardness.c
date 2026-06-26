/* Example: SETH Conditional Lower Bounds */
#include <stdio.h>
#include <stdlib.h>
#include "../include/seth_hardness.h"
#include "../include/frechet_distance.h"

int main(void) {
    printf("=== SETH Conditional Lower Bounds ===\n\n");

    printf("1. Lower bound status reports:\n");
    print_edit_seth_status();
    printf("\n");
    print_lcs_seth_status();
    printf("\n");
    print_frechet_seth_status();

    printf("\n2. Comprehensive report:\n");
    print_comprehensive_lower_bound_report();

    printf("\n3. OV problem demo:\n");
    ov_instance_t *inst = ov_instance_create(5, 5, 16);
    ov_instance_random_fill(inst, 42);
    ov_result_t res = ov_brute_force(inst);
    printf("   OV(5,5,16): orthogonal=%s, ops=%llu\n",
           res.exists_orthogonal ? "YES" : "NO",
           (unsigned long long)res.operations);
    ov_instance_destroy(inst);

    printf("\n4. SETH constant analysis:\n");
    for (int32_t k = 3; k <= 8; k++) {
        printf("   s_%d = %.4f\n", k, seth_constant_s_k(k));
    }

    printf("\n5. Verification of hardness theorems:\n");
    printf("   Backurs-Indyk (n=100, eps=0.1): %s\n",
           backurs_indyk_verify(100, 0.1) ? "REFUTES SETH" : "consistent");
    printf("   ABW (n=100, eps=0.1): %s\n",
           abw_verify(100, 0.1) ? "REFUTES SETH" : "consistent");
    printf("   Bringmann-Kunnemann (n=100, eps=0.1): %s\n",
           bringmann_kunnemann_verify(100, 0.1) ? "REFUTES SETH" : "consistent");

    printf("\n6. Frechet distance demo:\n");
    double px[] = {0, 1, 2, 3};
    double py[] = {0, 1, 0, 1};
    double qx[] = {0, 1, 2, 3};
    double qy[] = {0, 0.5, 0, 0.5};
    curve_t *c1 = curve_from_arrays(px, py, 4);
    curve_t *c2 = curve_from_arrays(qx, qy, 4);
    if (c1 && c2) {
        printf("   Frechet = %.4f\n", discrete_frechet_distance(c1, c2));
        printf("   DTW     = %.4f\n", dtw_distance(c1, c2));
        curve_destroy(c1); curve_destroy(c2);
    }

    printf("\n7. Breakthrough implications:\n");
    breakthrough_implications(1.5);
    printf("\n");
    breakthrough_implications(2.0);

    return 0;
}
