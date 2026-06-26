/* seth_empirical.c -- Empirical SETH Validation (L7, L8)
 *
 * Provides tools for empirically testing SETH predictions:
 * - Runtime extrapolation from small to large instances
 * - Checking whether observed exponents approach 1 as k grows
 * - Measuring the "effective" s_k for a given SAT solver
 * - Confidence intervals for SETH violation claims
 *
 * This connects the theoretical SETH/ETH framework to
 * actual computational experiments.
 */
#include "seth.h"
#include "cnf_generator.h"
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdio.h>
#include <math.h>

/* ===================================================================
 * Runtime Measurement and Extrapolation
 * =================================================================== */

/* Measure the effective exponent of a SAT solver on a given instance.
 * Runs the solver multiple times and computes:
 *   c = log2(average_ops) / n
 *
 * The effective exponent should be compared against seth_limit_s_k(k)
 * to check SETH consistency. */
typedef struct {
    double   mean_exponent;
    double   std_exponent;
    double   min_exponent;
    double   max_exponent;
    int32_t  num_trials;
    uint64_t mean_operations;
} runtime_stats_t;

runtime_stats_t measure_solver_exponent(
    sat_result_t (*solver)(const cnf_formula_t *),
    const cnf_formula_t *f,
    int32_t num_trials) {

    runtime_stats_t stats = {0.0, 0.0, INFINITY, -INFINITY, num_trials, 0};

    if (num_trials <= 0) num_trials = 5;
    if (num_trials > 100) num_trials = 100;

    double *exponents = (double *)malloc((size_t)num_trials * sizeof(double));
    if (!exponents) return stats;

    uint64_t total_ops = 0;

    for (int32_t t = 0; t < num_trials; t++) {
        sat_result_t result = solver(f);
        uint64_t ops = result.num_branches +
                       result.num_backtracks +
                       result.num_propagations;
        if (ops == 0) ops = 1; /* avoid log(0) */

        double exp = compute_sat_exponent(ops, f->num_vars);
        exponents[t] = exp;
        total_ops += ops;

        if (exp < stats.min_exponent) stats.min_exponent = exp;
        if (exp > stats.max_exponent) stats.max_exponent = exp;
    }

    /* Compute mean */
    double sum = 0.0;
    for (int32_t t = 0; t < num_trials; t++)
        sum += exponents[t];
    stats.mean_exponent = sum / (double)num_trials;

    /* Compute standard deviation */
    double sum_sq = 0.0;
    for (int32_t t = 0; t < num_trials; t++) {
        double diff = exponents[t] - stats.mean_exponent;
        sum_sq += diff * diff;
    }
    stats.std_exponent = sqrt(sum_sq / (double)num_trials);

    stats.mean_operations = total_ops / (uint64_t)num_trials;

    free(exponents);
    return stats;
}

/* ===================================================================
 * SETH Consistency Check
 *
 * Given measured exponents for increasing k, check if the trend
 * is consistent with SETH (s_k ? 1 as k ? ?).
 *
 * SETH predicts:
 *   s_3 ? 0.386
 *   s_4 ? 0.555
 *   s_5 ? 0.650
 *   ...
 *   s_? = 1.0
 *
 * If the measured s_k values plateau below 1.0 for all k,
 * that would be evidence AGAINST SETH (could be refuted
 * by finding an algorithm with exponent bounded away from 1).
 * =================================================================== */

typedef struct {
    int32_t  k;               /* clause width */
    double   measured_s_k;    /* effective exponent */
    double   theoretical_s_k; /* seth_limit_s_k(k) */
    double   gap;             /* theoretical - measured */
    bool     consistent;       /* measured >= theoretical (no violation) */
} seth_consistency_point_t;

/* Check SETH consistency across multiple k values.
 * For each k, measure the effective exponent of the solver
 * on random k-SAT instances and compare against seth_limit_s_k(k). */
void check_seth_consistency(
    sat_result_t (*solver)(const cnf_formula_t *),
    int32_t n,
    int32_t k_min,
    int32_t k_max,
    int32_t num_instances) {

    printf("\n=== SETH Consistency Check ===\n");
    printf("Solver on random k-SAT, n=%d, instances=%d\n", n, num_instances);
    printf("%4s %12s %12s %12s %10s\n",
           "k", "Measured s_k", "Theory s_k", "Gap", "Status");
    printf("----------------------------------------------------------\n");

    for (int32_t k = k_min; k <= k_max; k++) {
        double sum_exp = 0.0;
        int32_t valid = 0;

        for (int32_t inst = 0; inst < num_instances; inst++) {
            cnf_formula_t *f = generate_random_ksat(n, k,
                ksat_phase_transition(k));
            if (!f) continue;

            runtime_stats_t stats = measure_solver_exponent(solver, f, 3);
            sum_exp += stats.mean_exponent;
            valid++;

            cnf_destroy(f);
        }

        if (valid > 0) {
            double measured = sum_exp / (double)valid;
            double theory = seth_limit_s_k(k);
            double gap = theory - measured;
            bool consistent = measured >= theory - 0.05; /* tolerance */

            printf("%4d %12.4f %12.4f %12.4f %10s\n",
                   k, measured, theory, gap,
                   consistent ? "OK" : "WARNING");
        }
    }
    printf("----------------------------------------------------------\n");
    printf("WARNING = measured exponent below SETH prediction.\n");
    printf("This could indicate a SETH-violating algorithm.\n");
}

/* ===================================================================
 * Extrapolation: Predict Runtime on Larger Instances
 *
 * Given measured runtimes for small n, fit an exponential curve
 * T(n) = a * 2^{c*n} and predict runtime for larger n.
 *
 * This is how SETH-based lower bounds are empirically tested:
 * if T(n) grows like 2^{c*n} with c close to s_k, then SETH
 * predictions are consistent with the data.
 * =================================================================== */

/* Linear regression on log-transformed data:
 * log2(T) = log2(a) + c * n
 * Fit (n_i, log2(T_i)) for i=1..m to get c (slope). */

typedef struct {
    double   fitted_exponent;   /* slope c */
    double   intercept;         /* log2(a) */
    double   r_squared;         /* goodness of fit */
    double   predicted_time_100; /* predicted time for n=100 */
} extrapolation_result_t;

extrapolation_result_t extrapolate_runtime(
    sat_result_t (*solver)(const cnf_formula_t *),
    int32_t k,
    int32_t n_min,
    int32_t n_max,
    int32_t n_step) {

    extrapolation_result_t result = {0.0, 0.0, 0.0, 0.0};

    /* Collect data points */
    int32_t num_points = (n_max - n_min) / n_step + 1;
    if (num_points < 2) return result;

    double *n_vals = (double *)malloc((size_t)num_points * sizeof(double));
    double *log_t = (double *)malloc((size_t)num_points * sizeof(double));
    if (!n_vals || !log_t) {
        free(n_vals); free(log_t);
        return result;
    }

    int32_t point = 0;
    for (int32_t n = n_min; n <= n_max; n += n_step) {
        cnf_formula_t *f = generate_random_ksat(n, k,
            ksat_phase_transition(k));
        if (!f) continue;

        runtime_stats_t stats = measure_solver_exponent(solver, f, 3);

        n_vals[point] = (double)n;
        log_t[point] = log2((double)stats.mean_operations);
        point++;

        cnf_destroy(f);
    }

    if (point < 2) {
        free(n_vals); free(log_t);
        return result;
    }

    /* Linear regression: log_t = intercept + exponent * n */
    double sum_n = 0, sum_logt = 0, sum_n2 = 0, sum_n_logt = 0;
    for (int32_t i = 0; i < point; i++) {
        sum_n += n_vals[i];
        sum_logt += log_t[i];
        sum_n2 += n_vals[i] * n_vals[i];
        sum_n_logt += n_vals[i] * log_t[i];
    }

    double denom = (double)point * sum_n2 - sum_n * sum_n;
    if (fabs(denom) < 1e-10) {
        free(n_vals); free(log_t);
        return result;
    }

    result.intercept = (sum_n2 * sum_logt - sum_n * sum_n_logt) / denom;
    result.fitted_exponent = ((double)point * sum_n_logt - sum_n * sum_logt) / denom;

    /* R-squared */
    double mean_logt = sum_logt / (double)point;
    double ss_res = 0.0, ss_tot = 0.0;
    for (int32_t i = 0; i < point; i++) {
        double predicted = result.intercept + result.fitted_exponent * n_vals[i];
        ss_res += (log_t[i] - predicted) * (log_t[i] - predicted);
        ss_tot += (log_t[i] - mean_logt) * (log_t[i] - mean_logt);
    }
    result.r_squared = 1.0 - ss_res / ss_tot;

    /* Predict for n=100 */
    result.predicted_time_100 = pow(2.0,
        result.intercept + result.fitted_exponent * 100.0);

    free(n_vals);
    free(log_t);
    return result;
}

/* ===================================================================
 * SETH Evidence Assessment
 *
 * Given empirical data, assess the strength of evidence
 * for or against SETH.
 *
 * Evidence FOR SETH:
 * - Measured s_k increases monotonically with k
 * - s_k approaches 1 as k grows
 * - No algorithm achieves exponent < 1 - eps for all k
 *
 * Evidence AGAINST SETH:
 * - s_k plateaus below 1.0 for all measured k
 * - An algorithm achieves (1-eps) exponent for large k
 * - Breakthrough: new algorithm with exponent independent of k
 * =================================================================== */

typedef enum {
    SETH_SUPPORTED,       /* evidence consistent with SETH */
    SETH_WEAK_EVIDENCE,   /* insufficient data */
    SETH_TENSION,         /* some tension with data */
    SETH_REFUTED          /* strong evidence against SETH */
} seth_evidence_level_t;

seth_evidence_level_t assess_seth_evidence(
    const double *measured_s_k,
    int32_t k_min,
    int32_t k_max) {

    int32_t num_points = k_max - k_min + 1;
    if (num_points < 3) return SETH_WEAK_EVIDENCE;

    /* Check monotonicity: s_k should increase with k */
    bool monotonic = true;
    for (int32_t i = 1; i < num_points; i++) {
        if (measured_s_k[i] < measured_s_k[i-1] - 0.01)
            monotonic = false;
    }

    /* Check approach to 1: s_k should be close to 1 for large k */
    double last_s_k = measured_s_k[num_points - 1];
    bool approaching_one = (last_s_k > 0.9);

    /* Check for SETH violation: any s_k significantly below
     * the theoretical limit */
    bool any_violation = false;
    for (int32_t i = 0; i < num_points; i++) {
        int32_t k = k_min + i;
        double theory = seth_limit_s_k(k);
        if (measured_s_k[i] < theory - 0.1)
            any_violation = true;
    }

    if (any_violation)
        return SETH_REFUTED;
    else if (!approaching_one)
        return SETH_TENSION;
    else if (monotonic && approaching_one)
        return SETH_SUPPORTED;
    else
        return SETH_WEAK_EVIDENCE;
}

/* ===================================================================
 * Print SETH Status Report
 * =================================================================== */

void print_seth_status_report(void) {
    printf("\n");
    printf("==============================================\n");
    printf("     STRONG EXPONENTIAL TIME HYPOTHESIS       \n");
    printf("               Status Report                  \n");
    printf("==============================================\n\n");

    printf("SETH Statement:\n");
    printf("  For every epsilon > 0, there exists k such that\n");
    printf("  k-SAT cannot be solved in O(2^{(1-epsilon)n}) time.\n\n");

    printf("Current Best k-SAT Exponents (PPSZ bounds):\n");
    for (int32_t k = 3; k <= 10; k++) {
        double s = seth_limit_s_k(k);
        printf("  s_%d = %.4f  (base = 2^{%.4f n})\n", k, s, s);
    }
    printf("  s_inf = 1.0000  (SETH limit)\n\n");

    printf("ETH Status:\n");
    printf("  3-SAT exponent s_3 = %.4f\n", seth_limit_s_k(3));
    printf("  ETH holds if s_3 > 0 (currently %.4f > 0 ?)\n",
           seth_limit_s_k(3));

    printf("\nKey Consequences if SETH holds:\n");
    printf("  - OV requires n^{2-o(1)} time\n");
    printf("  - Edit Distance requires n^{2-o(1)} time\n");
    printf("  - LCS requires n^{2-o(1)} time\n");
    printf("  - k-Dominating Set requires n^{k-o(1)} time\n");
    printf("  - Frechet Distance requires n^{2-o(1)} time\n");

    printf("\nKey Consequences if SETH is false:\n");
    printf("  - P != NP still possible (ETH may still hold)\n");
    printf("  - Faster algorithms for all problems above\n");
    printf("  - CNF-SAT has 2^{(1-eps)n} algorithm for all eps\n");

    printf("\n==============================================\n");
}

/* ===================================================================
 * Runtime Prediction for Specific Problems under SETH
 * =================================================================== */

/* Given SETH, predict the minimum time needed for Edit Distance
 * on strings of length N.
 * Under SETH: T(N) = Omega(N^2).
 * Best known: O(N^2 / log^2 N) (Masek-Paterson 1980, Four Russians)
 *              O(N^2 / log^{100} N) improvements exist but no
 *              O(N^{2-eps}) algorithm. */
double seth_predict_edit_distance_time(int32_t N, double ops_per_second) {
    /* Under SETH: no O(N^{2-eps}) algorithm.
     * Lower bound: T(N) >= c * N^2 for some constant c. */
    double ops = (double)N * (double)N;
    return ops / ops_per_second;
}

/* Predict minimum time for Orthogonal Vectors.
 * Given N vectors of dimension d.
 * Under SETH: T(N, d) = N^{2-o(1)} for d = omega(log N).
 * Best known: O(N^2 / 2^{Omega(d)}) for small d (trivial improvement). */
double seth_predict_ov_time(int32_t N, int32_t d, double ops_per_second) {
    /* Under SETH: N^{2-o(1)} for d = omega(log N).
     * For d = c*log N: T = N^{2 - O(1/c)}. */
    (void)d;
    double ops = pow((double)N, 2.0);
    return ops / ops_per_second;
}

/* Predict minimum time for k-Dominating Set on n-vertex graph.
 * Under SETH: n^{k-o(1)}. */
double seth_predict_kdom_time(int32_t n, int32_t k, double ops_per_second) {
    double ops = pow((double)n, (double)k);
    return ops / ops_per_second;
}

/* Predict minimum time for APSP (All-Pairs Shortest Paths).
 * Under APSP conjecture: n^{3-o(1)}.
 * Best known: O(n^3 / 2^{Omega(sqrt(log n))}) (Williams 2014)
 *              or O(n^3 / log^{100} n) using combinatorial methods.
 * Matrix multiplication: O(n^{omega}) where omega < 2.373,
 * but this gives APSP only for small integer weights. */
double apsp_conjecture_time(int32_t n, double ops_per_second) {
    /* APSP conjecture: no O(n^{3-eps}) algorithm.
     * Lower bound: Omega(n^3). */
    double ops = pow((double)n, 3.0);
    return ops / ops_per_second;
}
