#ifndef COMPLEXITY_H
#define COMPLEXITY_H
#include <stdint.h>
#include <stdbool.h>
#include <time.h>
#include "threesum.h"
#include "apsp.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef struct { int32_t n; double time_seconds; int64_t operations; double ops_per_sec; } cx_measurement_t;
typedef struct { cx_measurement_t *points; int32_t count; double estimated_exponent; double r_squared; double constant_factor; } cx_scaling_result_t;
double cx_measure_time(void (*fn)(void*), void *arg, int32_t iterations);
cx_scaling_result_t *cx_estimate_exponent(cx_measurement_t *pts, int32_t count);
void cx_scaling_result_destroy(cx_scaling_result_t *r);
void cx_print_measurement(const cx_measurement_t *m);
void cx_print_scaling(const cx_scaling_result_t *r);
bool cx_is_subcubic(double exponent, int32_t n);
bool cx_is_subquadratic(double exponent, int32_t n);
void cx_bench_3sum_scaling(void);
void cx_bench_apsp_scaling(void);
void cx_print_all_benchmarks(void);
void cx_print_complexity_table(void);
void cx_self_test(void);
void cx_verify_conjecture_empirically(void);
void cx_run_full_analysis(void);
#ifdef __cplusplus
}
#endif
#endif

/*
 * Complexities tracked:
 *   O(n^3): APSP, Min-Plus Product (cubic barrier)
 *   O(n^2): 3SUM, OV, Edit Distance (quadratic barrier)
 *   O(2^n): SAT (exponential barrier - SETH)
 *
 * Key insight: Fine-grained complexity unifies these barriers
 * through a web of reductions showing that breaking one barrier
 * would break them all.
 */