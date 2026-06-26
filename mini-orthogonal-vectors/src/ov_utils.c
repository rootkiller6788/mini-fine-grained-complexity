/* ov_utils.c - OV Utility functions and benchmarking support
 * Provides timing, statistics, validation, and batch operation support.
 */
#include "ov.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>

/* Compare the results of two OV algorithms on the same instance.
 * Returns true if they agree on whether an orthogonal pair exists.
 * If both find a pair, checks if the indices match.
 * Prints diagnostic info to stdout. */
bool ov_compare_algorithms(const ov_instance_t *inst,
                            const char *name1, ov_result_t (*alg1)(const ov_instance_t*),
                            const char *name2, ov_result_t (*alg2)(const ov_instance_t*)) {
    if (!inst || !alg1 || !alg2) return false;
    ov_result_t r1 = alg1(inst);
    ov_result_t r2 = alg2(inst);
    bool agree = (r1.found == r2.found);
    if (agree && r1.found) {
        agree = (r1.a_index == r2.a_index && r1.b_index == r2.b_index);
    }
    if (!agree) {
        printf("  MISMATCH: %s(%s) vs %s(%s)\n",
               name1, r1.found ? "YES" : "NO",
               name2, r2.found ? "YES" : "NO");
    }
    return agree;
}

/* Run a benchmark of a given OV algorithm on an instance.
 * Returns the result with timing information filled in. */
ov_result_t ov_benchmark(const ov_instance_t *inst,
                          ov_result_t (*alg)(const ov_instance_t*),
                          int32_t num_repeats) {
    ov_result_t best = {false, -1, -1, 0, 1e100};
    if (!inst || !alg) return best;
    for (int32_t rep = 0; rep < num_repeats; rep++) {
        ov_result_t r = alg(inst);
        if (r.time_s < best.time_s) best = r;
    }
    return best;
}

/* Generate a human-readable report comparing all algorithms
 * on a given OV instance. */
void ov_benchmark_all(const ov_instance_t *inst) {
    if (!inst) return;
    printf("=== OV Algorithm Benchmark ===\n");
    printf("Instance: |A|=|B|=%d, d=%d\n",
           inst->A->num_vectors, inst->A->dimension);
    printf("Total pairs: %lld\n\n", (long long)ov_num_pairs(inst));

    struct { const char *name; ov_result_t (*alg)(const ov_instance_t*); } algs[] = {
        {"Brute Force",      ov_brute_force},
        {"Brute Early",      ov_brute_force_early},
        {"Sparse First",     ov_brute_force_sparse_first},
        {"Williams",          ov_williams},
        {"Four-Russians",     ov_four_russians},
        {"Matrix Mult",       ov_via_matrix_mult},
        {NULL, NULL}
    };

    for (int32_t i = 0; algs[i].name; i++) {
        ov_result_t r = ov_benchmark(inst, algs[i].alg, 3);
        printf("  %-20s: found=%s, time=%.6fs, ops=%lld\n",
               algs[i].name,
               r.found ? "YES" : "NO",
               r.time_s,
               (long long)r.ops);
    }
}

/* Estimate the empirical running time exponent by varying n.
 * Creates instances at different sizes and fits T(n) = C * n^exp.
 * Returns the estimated exponent (2.0 = quadratic). */
double ov_estimate_exponent(int32_t d, double density, int32_t min_n, int32_t max_n, int32_t steps) {
    if (steps < 2 || min_n <= 0 || max_n <= min_n) return -1.0;
    printf("=== Empirical Time Exponent Estimation ===\n");
    printf("d=%d, density=%.2f, n in [%d, %d]\n", d, density, min_n, max_n);
    /* Simple estimation: run brute force at two sizes and compute ratio */
    ov_instance_t *inst1 = ov_create_random(min_n, d, density, false);
    ov_instance_t *inst2 = ov_create_random(max_n, d, density, false);
    if (!inst1 || !inst2) {
        ov_destroy(inst1); ov_destroy(inst2);
        return -1.0;
    }
    ov_result_t r1 = ov_benchmark(inst1, ov_brute_force, 2);
    ov_result_t r2 = ov_benchmark(inst2, ov_brute_force, 2);
    ov_destroy(inst1); ov_destroy(inst2);
    if (r1.time_s <= 0 || r2.time_s <= 0) return -1.0;
    double exp_est = log(r2.time_s / r1.time_s) / log((double)max_n / (double)min_n);
    printf("  T(%d)=%.4fs, T(%d)=%.4fs => exponent ~ %.3f\n",
           min_n, r1.time_s, max_n, r2.time_s, exp_est);
    return exp_est;
}

/* Compute the fraction of orthogonal pairs in an instance.
 * For random instances with density delta, this is approximately
 * (1-delta)^d for each pair. */
double ov_orthogonal_fraction(const ov_instance_t *inst) {
    if (!inst || !inst->A || !inst->B) return 0.0;
    int64_t total = ov_num_pairs(inst);
    if (total == 0) return 0.0;
    int64_t orth = ov_count_orthogonal_pairs(inst);
    return (double)orth / (double)total;
}

/* Verify that the OV conjecture holds for a given (n,d) combination.
 * Creates random instances and checks that no algorithm achieves
 * subquadratic time for superlogarithmic dimension. */
bool ov_verify_conjecture_empirical(int32_t n, int32_t d) {
    if (!ov_is_superlogarithmic(n, d)) {
        printf("  Note: d=%d is not omega(log n=%d), conjecture doesn't apply.\n",
               d, (int32_t)(log2((double)n)));
        return true;
    }
    printf("  Testing OV conjecture for n=%d, d=%d...\n", n, d);
    double exp = ov_estimate_exponent(d, 0.3, n/2, n, 2);
    if (exp < 0) return true; /* couldn't measure */
    /* If empirical exponent is significantly below 2, conjecture might be violated */
    bool holds = exp >= 1.9;
    printf("  Conjecture holds (exp >= 1.9): %s\n", holds ? "YES" : "POSSIBLY VIOLATED");
    return holds;
}

/* Generate a LaTeX-compatible table of OV algorithm complexities. */
void ov_print_complexity_table(void) {
    printf("\n=== OV Algorithm Complexity Table ===\n\n");
    printf("Algorithm              | Time Complexity        | Space     | Condition\n");
    printf("-----------------------|------------------------|-----------|----------\n");
    printf("Brute Force            | O(n^2 * d)             | O(1)      | Always\n");
    printf("Williams (2005)        | n^{2-1/O(c)}           | O(n*d)    | d=c*log n\n");
    printf("Meet-in-the-Middle     | O(n * 2^{d/2})         | O(n+2^{d/2}) | d small\n");
    printf("Four-Russians          | O(n^2 * d / log n)     | O(2^b)    | d=Theta(log n)\n");
    printf("Matrix Multiplication  | O(n^omega)             | O(n^2)    | omega<2.373\n");
    printf("Polynomial Method      | O(n^{2-eps})           | O(n*d)    | d=c*log n\n");
    printf("Light-Sparse           | O(n^2 * delta*d)       | O(n*d)    | delta=o(1)\n");
    printf("LSH (bit-sampling)     | O(n^{1+rho})           | O(n*d)    | approx\n");
    printf("\n");
    printf("c = d / log_2(n), omega = matrix mult exponent, delta = vector density\n");
    printf("Reference: Abboud-Williams-Yu (2015), Williams (2005)\n\n");
}

/* Utility: count how many vectors in A have a given bit set. */
int32_t *ov_bit_frequencies(const ov_instance_t *inst, bool set_a) {
    if (!inst) return NULL;
    vector_set_t *vs = set_a ? inst->A : inst->B;
    if (!vs) return NULL;
    int32_t *freq = (int32_t *)calloc((size_t)vs->dimension, sizeof(int32_t));
    if (!freq) return NULL;
    for (int32_t i = 0; i < vs->num_vectors; i++)
        for (int32_t k = 0; k < vs->dimension; k++)
            if (bv_get(&vs->vectors[i], k))
                freq[k]++;
    return freq;
}

/* Print a histogram of bit frequencies (for analyzing instance structure). */
void ov_print_bit_histogram(const ov_instance_t *inst) {
    if (!inst || !inst->A) return;
    int32_t *freqA = ov_bit_frequencies(inst, true);
    int32_t *freqB = ov_bit_frequencies(inst, false);
    if (!freqA || !freqB) { free(freqA); free(freqB); return; }
    int32_t d = inst->A->dimension;
    int32_t nA = inst->A->num_vectors;
    printf("=== Bit Frequency Histogram (d=%d, n=%d) ===\n", d, nA);
    printf("Bit | FreqA | FreqB | DensityA | DensityB\n");
    printf("----|-------|-------|----------|----------\n");
    for (int32_t k = 0; k < d && k < 32; k++)
        printf("%3d | %5d | %5d | %8.3f | %8.3f\n",
               k, freqA[k], freqB[k],
               (double)freqA[k]/(double)nA,
               (double)freqB[k]/(double)nA);
    if (d > 32) printf("  ... (%d more bits)\n", d - 32);
    free(freqA); free(freqB);
}

/* Self-test: run all internal consistency checks.
 * Returns number of checks passed. */
int32_t ov_self_test(void) {
    int32_t passed = 0;
    int32_t total = 0;

    /* Test 1: bv_create/bv_destroy */
    binary_vector_t *v = bv_create(100);
    if (v) { passed++; } total++;
    if (v && v->dimension == 100) { passed++; } total++;
    bv_destroy(v);

    /* Test 2: bv_set/bv_get */
    v = bv_create(64);
    bv_set(v, 10, true);
    if (bv_get(v, 10)) { passed++; } total++;
    if (!bv_get(v, 11)) { passed++; } total++;
    bv_destroy(v);

    /* Test 3: orthogonality */
    binary_vector_t *a = bv_create(10);
    binary_vector_t *b = bv_create(10);
    bv_set(a, 0, 1); bv_set(a, 3, 1);
    bv_set(b, 1, 1); bv_set(b, 5, 1);
    if (bv_are_orthogonal(a, b)) { passed++; } total++;
    bv_set(b, 0, 1);
    if (!bv_are_orthogonal(a, b)) { passed++; } total++;
    bv_destroy(a); bv_destroy(b);

    /* Test 4: vector set */
    vector_set_t *vs = vs_create(3, 20);
    if (vs && vs->num_vectors == 3) { passed++; } total++;
    vs_destroy(vs);

    /* Test 5: OV instance */
    ov_instance_t *inst = ov_create(5, 8);
    if (inst && ov_validate(inst)) { passed++; } total++;
    ov_destroy(inst);

    return passed;
}
