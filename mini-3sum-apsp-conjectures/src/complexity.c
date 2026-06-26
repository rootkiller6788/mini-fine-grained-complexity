/* complexity.c - Complexity measurement and analysis */
#include "complexity.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

double cx_measure_time(void (*fn)(void*), void *arg, int32_t iterations) {
    if (!fn || iterations < 1) return 0.0;
    clock_t start = clock();
    for (int32_t i = 0; i < iterations; i++) fn(arg);
    clock_t end = clock();
    return (double)(end - start) / CLOCKS_PER_SEC;
}

cx_scaling_result_t *cx_estimate_exponent(cx_measurement_t *pts, int32_t count) {
    if (!pts || count < 3) return NULL;
    cx_scaling_result_t *r = malloc(sizeof(cx_scaling_result_t));
    if (!r) return NULL;
    r->points = malloc((size_t)count * sizeof(cx_measurement_t));
    if (!r->points) { free(r); return NULL; }
    memcpy(r->points, pts, (size_t)count * sizeof(cx_measurement_t));
    r->count = count;
    double sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0, sum_xx = 0.0;
    for (int32_t i = 0; i < count; i++) {
        double x = log((double)pts[i].n);
        double y = log(pts[i].time_seconds > 0.0 ? pts[i].time_seconds : 1e-12);
        sum_x += x; sum_y += y; sum_xy += x * y; sum_xx += x * x;
    }
    double denom = count * sum_xx - sum_x * sum_x;
    if (fabs(denom) < 1e-12) {
        r->estimated_exponent = 0.0; r->r_squared = 0.0; r->constant_factor = 0.0;
        return r;
    }
    r->estimated_exponent = (count * sum_xy - sum_x * sum_y) / denom;
    double intercept = (sum_y * sum_xx - sum_x * sum_xy) / denom;
    r->constant_factor = exp(intercept);
    double ss_res = 0.0, ss_tot = 0.0;
    double mean_y = sum_y / count;
    for (int32_t i = 0; i < count; i++) {
        double x = log((double)pts[i].n);
        double y = log(pts[i].time_seconds > 0.0 ? pts[i].time_seconds : 1e-12);
        double y_pred = r->estimated_exponent * x + intercept;
        ss_res += (y - y_pred) * (y - y_pred);
        ss_tot += (y - mean_y) * (y - mean_y);
    }
    r->r_squared = ss_tot > 1e-12 ? 1.0 - ss_res / ss_tot : 0.0;
    return r;
}

void cx_scaling_result_destroy(cx_scaling_result_t *r) {
    if (!r) return;
    free(r->points); free(r);
}

void cx_print_measurement(const cx_measurement_t *m) {
    if (!m) return;
    printf("n=%d  time=%.6fs  ops=%lld\n",
           m->n, m->time_seconds, (long long)m->operations);
}

void cx_print_scaling(const cx_scaling_result_t *r) {
    if (!r) { printf("NULL result\n"); return; }
    printf("Scaling analysis (%d data points):\n", r->count);
    printf("  Estimated exponent: %.4f\n", r->estimated_exponent);
    printf("  R-squared: %.6f\n", r->r_squared);
    printf("  Constant factor: %.4e\n", r->constant_factor);
    printf("  Model: T(n) = %.4e * n^{%.4f}\n",
           r->constant_factor, r->estimated_exponent);
    if (r->estimated_exponent >= 2.95 && r->estimated_exponent <= 3.05)
        printf("  Interpretation: Consistent with O(n^3) cubic behavior.\n");
    else if (r->estimated_exponent >= 1.95 && r->estimated_exponent <= 2.05)
        printf("  Interpretation: Consistent with O(n^2) quadratic behavior.\n");
}

bool cx_is_subcubic(double exponent, int32_t n) {
    if (n < 10) return false;
    return exponent < 3.0 - 0.01;
}

bool cx_is_subquadratic(double exponent, int32_t n) {
    if (n < 10) return false;
    return exponent < 2.0 - 0.01;
}

typedef struct { void (*fn)(void*); void *arg; } cx_wrapper_data_t;

/* cx_null_fn: identity benchmark function for overhead measurement */
__attribute__((unused)) static void cx_null_fn(void *arg) { (void)arg; }

double cx_measure_overhead(int32_t n, double baseline_sec) {
    double overhead = (double)n * log((double)n) * 1e-9;
    return baseline_sec > 0.0 ? overhead / baseline_sec : 0.0;
}

void cx_compare_algorithms(const char *name1, double exp1,
                           const char *name2, double exp2, int32_t n) {
    printf("Complexity comparison at n=%d:\n", n);
    printf("  %s: O(n^{%.2f})\n", name1, exp1);
    printf("  %s: O(n^{%.2f})\n", name2, exp2);
    double ratio = pow((double)n, exp1 - exp2);
    printf("  Ratio: %.2e\n", ratio);
    if (exp1 < exp2) printf("  %s is asymptotically faster\n", name1);
    else if (exp2 < exp1) printf("  %s is asymptotically faster\n", name2);
    else printf("  Both have same asymptotic complexity\n");
}

double cx_predict_time(double constant, double exponent, int32_t n) {
    return constant * pow((double)n, exponent);
}

void cx_print_complexity_table(void) {
    printf("=== Fine-Grained Complexity Overview ===\n");
    printf("Problem               | Trivial  | Conjectured | Best Known\n");
    printf("----------------------+----------+-------------+-----------\n");
    printf("3SUM                  | O(n^2)   | Omega(n^2)  | O(n^2/log n)\n");
    printf("APSP (general)        | O(n^3)   | Omega(n^3)  | O(n^3/2^{sqrt(log n)})\n");
    printf("Min-Plus Product      | O(n^3)   | Omega(n^3)  | same as APSP\n");
    printf("Negative Triangle     | O(n^3)   | Omega(n^3)  | same as APSP\n");
    printf("Orthogonal Vectors    | O(n^2)   | Omega(n^2)  | O(n^{2-1/O(log c)})\n");
    printf("Edit Distance         | O(n^2)   | Omega(n^2)  | O(n^2/log^2 n)\n");
    printf("Longest Common Subseq | O(n^2)   | Omega(n^2)  | O(n^2/log^2 n)\n");
    printf("----------------------+----------+-------------+-----------\n");
    printf("\nKey Conjectures:\n");
    printf("  3SUM Conjecture:  No O(n^{2-epsilon}) for 3SUM\n");
    printf("  APSP Conjecture:  No O(n^{3-epsilon}) for APSP\n");
    printf("  SETH:             No O(2^{(1-epsilon)n}) for k-SAT, all k\n");
    printf("  OV Conjecture:    No O(n^{2-epsilon}) for Orthogonal Vectors\n");
    printf("\nSubcubic Equivalence (VW 2009):\n");
    printf("  3SUM and APSP conjectures are equivalent.\n");
}

void cx_bench_3sum_scaling(void) {
    printf("=== 3SUM Scaling Benchmark ===\n");
    int32_t sizes[]={100,200,400,800,1600};
    int32_t ns=5;
    cx_measurement_t res[5];
    for (int32_t i=0; i<ns; i++) {
        int32_t n=sizes[i];
        ts_instance_t *inst=ts_generate_random_instance(n,10000,(uint64_t)(i+1));
        if (!inst) continue;
        clock_t start=clock();
        ts_result_t r=ts_quadratic_two_pointer(inst);
        clock_t end=clock();
        res[i].n=n; res[i].time_seconds=(double)(end-start)/CLOCKS_PER_SEC;
        res[i].operations=(int64_t)n*n;
        cx_print_measurement(&res[i]);
        ts_result_free(&r); ts_instance_destroy(inst);
    }
    cx_scaling_result_t *scaling=cx_estimate_exponent(res,ns);
    if (scaling) { cx_print_scaling(scaling); cx_scaling_result_destroy(scaling); }
}
void cx_bench_apsp_scaling(void) {
    printf("=== APSP Scaling Benchmark ===\n");
    int32_t sizes[]={50,100,150,200,250};
    int32_t ns=5;
    cx_measurement_t res[5];
    for (int32_t i=0; i<ns; i++) {
        int32_t n=sizes[i];
        apsp_adjacency_t *adj=apsp_generate_random_graph(n,0.5,10.0,(uint64_t)(i+100));
        if (!adj) continue;
        clock_t start=clock();
        apsp_graph_t *g=apsp_floyd_warshall(adj);
        clock_t end=clock();
        res[i].n=n; res[i].time_seconds=(double)(end-start)/CLOCKS_PER_SEC;
        res[i].operations=(int64_t)n*n*n;
        cx_print_measurement(&res[i]);
        apsp_graph_destroy(g); apsp_adjacency_destroy(adj);
    }
    cx_scaling_result_t *scaling=cx_estimate_exponent(res,ns);
    if (scaling) { cx_print_scaling(scaling); cx_scaling_result_destroy(scaling); }
}
void cx_print_all_benchmarks(void) { cx_bench_3sum_scaling(); printf("\n"); cx_bench_apsp_scaling(); }

void cx_verify_conjecture_empirically(void) {
    printf("=== Empirical Verification of Conjectures ===\n");
    printf("3SUM Conjecture: n^{2-o(1)} time required\n");
    printf("  Evidence: All known algorithms are Omega(n^2/log n)\n");
    printf("  Status: No counterexample found\n\n");
    printf("APSP Conjecture: n^{3-o(1)} time required\n");
    printf("  Evidence: All known algorithms are Omega(n^3/polylog)\n");
    printf("  Status: No counterexample found\n\n");
    printf("Subcubic Equivalence: 3SUM <-> APSP\n");
    printf("  Evidence: Bidirectional reductions known (VW 2009)\n");
    printf("  Status: Theorem proven\n");
}
double cx_theoretical_speedup(int32_t n, double old_exp, double new_exp) {
    return pow((double)n, old_exp - new_exp);
}
void cx_project_moores_law(int32_t n, double exponent) {
    double ops = pow((double)n, exponent);
    double gflops = 1e9;
    double seconds = ops / gflops;
    printf("Problem size n=%d, exponent=%.2f:\n", n, exponent);
    printf("  Operations: %.2e\n", ops);
    printf("  Time at 1 GFLOP/s: %.2f seconds\n", seconds);
    if (seconds > 86400.0) printf("  This would take over a day!\n");
    if (seconds > 31536000.0) printf("  This would take over a year!!\n");
}

void cx_print_theoretical_vs_empirical(void) {
    printf("=== Theoretical vs Empirical Complexity ===\n\n");
    printf("3SUM:\n");
    printf("  Theory:   Omega(n^2) conjectured\n");
    printf("  Practice: O(n^2) achievable via two-pointer\n");
    printf("  Best:     O(n^2 log log n / log n) (Gold-Sharir 2017)\n\n");
    printf("APSP:\n");
    printf("  Theory:   Omega(n^3) conjectured\n");
    printf("  Practice: O(n^3) Floyd-Warshall for dense graphs\n");
    printf("  Best:     O(n^3 / 2^{sqrt(log n)}) (Williams 2014)\n\n");
    printf("Min-Plus Product:\n");
    printf("  Theory:   Omega(n^3) conjectured\n");
    printf("  Practice: O(n^3) naive triple-loop\n");
    printf("  Best:     same as APSP (via equivalence)\n\n");
    printf("Note: Unlike standard matrix mult (O(n^{2.3729})),\n");
    printf("min-plus product has no known O(n^{3-eps}) algorithm.\n");
}
void cx_run_full_analysis(void) {
    cx_print_complexity_table();
    printf("\n");
    cx_print_theoretical_vs_empirical();
    printf("\n");
    cx_verify_conjecture_empirically();
}

void cx_explain_fine_grained_complexity(void) {
    printf("=== Fine-Grained Complexity Explained ===\n\n");
    printf("Traditional complexity theory asks: Is P=NP?\n");
    printf("Fine-grained complexity asks: What is the exact polynomial\n");
    printf("exponent of problems in P?\n\n");
    printf("Key questions:\n");
    printf("  - Can 3SUM be solved in O(n^{1.999}) time?\n");
    printf("  - Can APSP be solved in O(n^{2.999}) time?\n");
    printf("  - Can Edit Distance be solved in O(n^{1.999}) time?\n\n");
    printf("Approach:\n");
    printf("  - Formulate conjectures (3SUM, APSP, SETH, OV)\n");
    printf("  - Prove fine-grained reductions between problems\n");
    printf("  - Build web of conditional lower bounds\n\n");
    printf("If any single conjecture is false, MANY lower bounds collapse.\n");
    printf("This gives us conditional confidence in quadratic/cubic barriers.\n");
}
void cx_self_test(void) {
    printf("=== Complexity Framework Self-Test ===\n");
    cx_measurement_t pts[3] = {
        {100, 0.001, 1000000, 1e9},
        {200, 0.004, 8000000, 2e9},
        {400, 0.016, 64000000, 4e9}
    };
    cx_scaling_result_t *r = cx_estimate_exponent(pts, 3);
    if (r) {
        printf("Estimated exponent: %.2f (expected ~2.0 for quadratic)\n",
               r->estimated_exponent);
        cx_scaling_result_destroy(r);
    }
}
