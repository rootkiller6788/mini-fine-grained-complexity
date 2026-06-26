/* ov_analysis.c - OV Analysis, Statistics and Profiling Tools
 * Provides tools for analyzing OV instance structure, algorithm profiling,
 * and empirical verification of theoretical predictions.
 */
#include "ov.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ===================================================================
 * Algorithm Profiling
 * =================================================================== */

/* Profile an OV algorithm: run it multiple times and collect statistics.
 * Returns the median running time (more robust than mean). */
typedef struct {
    double min_time, max_time, mean_time, median_time, stddev_time;
    int32_t num_runs;
    bool    all_agree;      /* whether all runs returned same result */
    int64_t total_ops;
} ov_profile_t;

ov_profile_t ov_profile_algorithm(const ov_instance_t *inst,
                                   ov_result_t (*alg)(const ov_instance_t*),
                                   int32_t num_runs) {
    ov_profile_t p;
    memset(&p, 0, sizeof(p));
    p.num_runs = num_runs;
    if (!inst || !alg || num_runs <= 0) return p;
    double *times = (double *)malloc((size_t)num_runs * sizeof(double));
    if (!times) return p;
    ov_result_t first = alg(inst);
    times[0] = first.time_s;
    p.total_ops = first.ops;
    p.min_time = times[0];
    p.max_time = times[0];
    double sum = times[0];
    p.all_agree = true;
    for (int32_t r = 1; r < num_runs; r++) {
        ov_result_t res = alg(inst);
        times[r] = res.time_s;
        sum += times[r];
        p.total_ops += res.ops;
        if (times[r] < p.min_time) p.min_time = times[r];
        if (times[r] > p.max_time) p.max_time = times[r];
        if (res.found != first.found) p.all_agree = false;
    }
    p.mean_time = sum / (double)num_runs;
    /* Compute median */
    for (int32_t i = 0; i < num_runs - 1; i++)
        for (int32_t j = i + 1; j < num_runs; j++)
            if (times[i] > times[j]) {
                double t = times[i]; times[i] = times[j]; times[j] = t;
            }
    if (num_runs % 2 == 0)
        p.median_time = (times[num_runs/2 - 1] + times[num_runs/2]) / 2.0;
    else
        p.median_time = times[num_runs/2];
    /* Compute stddev */
    double sum_sq = 0.0;
    for (int32_t r = 0; r < num_runs; r++) {
        double diff = times[r] - p.mean_time;
        sum_sq += diff * diff;
    }
    p.stddev_time = sqrt(sum_sq / (double)num_runs);
    free(times);
    return p;
}

void ov_print_profile(const char *algo_name, const ov_profile_t *p) {
    if (!algo_name || !p) return;
    printf("=== Profile: %s ===\n", algo_name);
    printf("  Runs: %d\n", p->num_runs);
    printf("  Min time:    %.6f s\n", p->min_time);
    printf("  Max time:    %.6f s\n", p->max_time);
    printf("  Mean time:   %.6f s\n", p->mean_time);
    printf("  Median time: %.6f s\n", p->median_time);
    printf("  Stddev:      %.6f s\n", p->stddev_time);
    printf("  All agree:   %s\n", p->all_agree ? "YES" : "NO");
    printf("  Total ops:   %lld\n\n", (long long)p->total_ops);
}

/* ===================================================================
 * Instance Structure Analysis
 * =================================================================== */

/* Compute the orthogonality matrix (full n x n matrix of dot products).
 * Returns the matrix as a flattened int32_t array of size n*n.
 * Caller must free. */
int32_t *ov_compute_dot_product_matrix(const ov_instance_t *inst) {
    if (!inst || !inst->A || !inst->B) return NULL;
    int32_t n = inst->A->num_vectors;
    int32_t *mat = (int32_t *)malloc((size_t)n * (size_t)n * sizeof(int32_t));
    if (!mat) return NULL;
    for (int32_t i = 0; i < n; i++)
        for (int32_t j = 0; j < n; j++)
            mat[i * n + j] = bv_dot_product(&inst->A->vectors[i],
                                             &inst->B->vectors[j]);
    return mat;
}

/* Compute statistics about the dot product distribution.
 * Useful for understanding instance structure beyond binary YES/NO. */
typedef struct {
    double  mean_dot;
    double  stddev_dot;
    int32_t min_dot;
    int32_t max_dot;
    int64_t num_zero_dots;    /* count of orthogonal pairs */
    double  zero_fraction;
} ov_dot_stats_t;

ov_dot_stats_t ov_compute_dot_statistics(const ov_instance_t *inst) {
    ov_dot_stats_t s;
    memset(&s, 0, sizeof(s));
    if (!inst || !inst->A || !inst->B) return s;
    int32_t n = inst->A->num_vectors;
    int32_t d = inst->A->dimension;
    int64_t total = (int64_t)n * (int64_t)n;
    s.min_dot = d + 1;
    s.max_dot = -1;
    double sum = 0.0;
    for (int32_t i = 0; i < n; i++) {
        for (int32_t j = 0; j < n; j++) {
            int32_t dot = bv_dot_product(&inst->A->vectors[i],
                                          &inst->B->vectors[j]);
            if (dot < s.min_dot) s.min_dot = dot;
            if (dot > s.max_dot) s.max_dot = dot;
            if (dot == 0) s.num_zero_dots++;
            sum += (double)dot;
        }
    }
    s.mean_dot = sum / (double)total;
    s.zero_fraction = (double)s.num_zero_dots / (double)total;
    /* Compute stddev */
    double sum_sq = 0.0;
    for (int32_t i = 0; i < n; i++)
        for (int32_t j = 0; j < n; j++) {
            double diff = (double)bv_dot_product(&inst->A->vectors[i],
                                                  &inst->B->vectors[j]) - s.mean_dot;
            sum_sq += diff * diff;
        }
    s.stddev_dot = sqrt(sum_sq / (double)total);
    return s;
}

void ov_print_dot_statistics(const ov_dot_stats_t *s) {
    if (!s) return;
    printf("=== Dot Product Distribution ===\n");
    printf("  Mean:    %.3f\n", s->mean_dot);
    printf("  Stddev:  %.3f\n", s->stddev_dot);
    printf("  Min:     %d\n", s->min_dot);
    printf("  Max:     %d\n", s->max_dot);
    printf("  #Zero:   %lld\n", (long long)s->num_zero_dots);
    printf("  %%Zero:   %.4f%%\n\n", s->zero_fraction * 100.0);
}

/* ===================================================================
 * Empirical SETH Verification
 * =================================================================== */

/* Run a parameter sweep and check if OV conjecture predictions hold.
 * For each (n,d) pair with d = omega(log n), verify that brute force
 * time scales quadratically in n. */
void ov_empirical_seth_check(int32_t d, double density,
                              int32_t n_start, int32_t n_end, int32_t n_step) {
    printf("\n=== Empirical SETH/OV Check ===\n");
    printf("d=%d, density=%.2f, n in [%d,%d] step %d\n\n",
           d, density, n_start, n_end, n_step);
    printf("n     | log n | d/log n | regime | brute_t     | predicted\n");
    printf("------|-------|---------|--------|-------------|----------\n");
    double prev_time = -1.0;
    int32_t prev_n = 0;
    for (int32_t n = n_start; n <= n_end; n += n_step) {
        ov_instance_t *inst = ov_create_random(n, d, density, false);
        if (!inst) continue;
        ov_result_t r = ov_benchmark(inst, ov_brute_force, 2);
        double log_n = log2((double)n);
        double ratio = (double)d / log_n;
        int32_t regime = ov_regime_classify(n, d);
        double exp_pred = ov_lower_bound_exponent(n, d);
        if (prev_time > 0 && prev_n > 0) {
            double emp_exp = log(r.time_s / prev_time) / log((double)n / (double)prev_n);
            printf("%-5d | %-5.1f | %-7.2f | %-6d | %-11.6f | exp=%.3f\n",
                   n, log_n, ratio, regime, r.time_s, emp_exp);
        } else {
            printf("%-5d | %-5.1f | %-7.2f | %-6d | %-11.6f | pred=%.3f\n",
                   n, log_n, ratio, regime, r.time_s, exp_pred);
        }
        prev_time = r.time_s;
        prev_n = n;
        ov_destroy(inst);
    }
    printf("\n");
}

/* ===================================================================
 * Visual ASCII Art for OV Instance
 * =================================================================== */

/* Print a small OV instance as an ASCII matrix visualization.
 * Only practical for small n,d (<= 20 each). */
void ov_print_ascii_visualization(const ov_instance_t *inst) {
    if (!inst || !inst->A || !inst->B) return;
    int32_t n = inst->A->num_vectors;
    int32_t d = inst->A->dimension;
    if (n > 20 || d > 20) {
        printf("Visualization skipped: n=%d, d=%d too large (>20)\n", n, d);
        return;
    }
    printf("\n=== OV Instance Visualization ===\n");
    printf("Set A (%d vectors x %d dims):\n", n, d);
    for (int32_t i = 0; i < n; i++) {
        printf("  A[%2d]: ", i);
        for (int32_t k = 0; k < d; k++)
            putchar(bv_get(&inst->A->vectors[i], k) ? '1' : '.');
        printf("  (wt=%d)\n", bv_hamming_weight(&inst->A->vectors[i]));
    }
    printf("Set B (%d vectors x %d dims):\n", n, d);
    for (int32_t j = 0; j < n; j++) {
        printf("  B[%2d]: ", j);
        for (int32_t k = 0; k < d; k++)
            putchar(bv_get(&inst->B->vectors[j], k) ? '1' : '.');
        printf("  (wt=%d)\n", bv_hamming_weight(&inst->B->vectors[j]));
    }
    /* Orthogonality matrix */
    printf("\nOrthogonality Matrix (n x n):\n    ");
    for (int32_t j = 0; j < n; j++) printf("%2d", j);
    printf("\n");
    for (int32_t i = 0; i < n; i++) {
        printf(" %2d ", i);
        for (int32_t j = 0; j < n; j++) {
            bool orth = bv_are_orthogonal(&inst->A->vectors[i],
                                           &inst->B->vectors[j]);
            printf(" %c", orth ? 'O' : '.');
        }
        printf("\n");
    }
    printf("  O = orthogonal pair found\n\n");
}

/* ===================================================================
 * Theoretical Bounds Calculator
 * =================================================================== */

/* Compute the theoretical maximum number of orthogonal pairs possible
 * for given n and d. This is an upper bound useful for sanity checks. */
int64_t ov_max_possible_orthogonal_pairs(int32_t n, int32_t d) {
    /* Upper bound: all n*n pairs could be orthogonal if vectors
     * are constructed with disjoint supports. But with d dimensions,
     * at most 2^d distinct vectors exist. For large n relative to d,
     * the pigeonhole principle limits orthogonal pairs.
     * Realistic bound: n^2 for constructive instances. */
    (void)d;
    return (int64_t)n * (int64_t)n;
}

/* Compute the minimum guaranteed number of orthogonal pairs
 * for any instance with given parameters.
 * This is 0 in general (can always construct a no-instance). */
int64_t ov_min_guaranteed_orthogonal_pairs(int32_t n, int32_t d) {
    (void)n; (void)d;
    return 0;  /* Always possible to have zero orthogonal pairs */
}

/* Given n and d, estimate the probability that a random OV instance
 * (with density 0.5) has at least one orthogonal pair.
 * This is the phase transition probability. */
double ov_phase_transition_probability(int32_t n, int32_t d) {
    /* Expected number of orthogonal pairs: n^2 * (3/4)^d
     * (since Pr[<a,b>=0] = (3/4)^d for density=0.5 independent bits).
     * By Poisson approximation: Pr[>=1 pair] ~ 1 - exp(-lambda)
     * where lambda = n^2 * (0.75)^d */
    double lambda = (double)n * (double)n * pow(0.75, (double)d);
    return 1.0 - exp(-lambda);
}

void ov_print_phase_transition_table(int32_t n_start, int32_t n_end, int32_t n_step,
                                       int32_t d_start, int32_t d_end, int32_t d_step) {
    printf("\n=== OV Phase Transition Table ===\n");
    printf("Probability that a random instance has >=1 orthogonal pair.\n\n");
    printf("n\d  ");
    for (int32_t d = d_start; d <= d_end; d += d_step)
        printf("%-8d", d);
    printf("\n");
    for (int32_t n = n_start; n <= n_end; n += n_step) {
        printf("%-4d ", n);
        for (int32_t d = d_start; d <= d_end; d += d_step) {
            double p = ov_phase_transition_probability(n, d);
            printf("%-8.4f", p);
        }
        printf("\n");
    }
    printf("\nDensity = 0.5 (worst-case for orthogonality probability)\n\n");
}
