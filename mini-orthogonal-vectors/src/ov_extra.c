/* ov_extra.c - Extended OV operations, GF(2) algebra, batch processing */
#include "ov.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

void bv_add_gf2(const binary_vector_t *a, const binary_vector_t *b, binary_vector_t *r) {
    bv_xor(a, b, r);
}

void bv_scalar_mul_gf2(const binary_vector_t *v, bool s, binary_vector_t *r) {
    if (!v || !r || v->dimension != r->dimension) return;
    if (s) bv_copy_into(v, r);
    else   bv_clear(r);
}

void bv_linear_combination_gf2(const binary_vector_t *a, bool c1,
                                const binary_vector_t *b, bool c2,
                                binary_vector_t *r) {
    if (!a || !b || !r) return;
    binary_vector_t *t1 = bv_create(a->dimension);
    binary_vector_t *t2 = bv_create(a->dimension);
    if (!t1 || !t2) { bv_destroy(t1); bv_destroy(t2); return; }
    bv_scalar_mul_gf2(a, c1, t1);
    bv_scalar_mul_gf2(b, c2, t2);
    bv_add_gf2(t1, t2, r);
    bv_destroy(t1); bv_destroy(t2);
}

bool bv_linearly_independent_gf2(const binary_vector_t *a, const binary_vector_t *b) {
    if (!a || !b) return false;
    if (bv_is_zero(a) || bv_is_zero(b)) return false;
    return !bv_equals(a, b);
}

int32_t vs_rank_gf2(const vector_set_t *vs) {
    if (!vs || vs->num_vectors == 0) return 0;
    int32_t n = vs->num_vectors, d = vs->dimension;
    binary_vector_t **mat = (binary_vector_t **)malloc((size_t)n * sizeof(binary_vector_t *));
    if (!mat) return 0;
    for (int32_t i = 0; i < n; i++) {
        mat[i] = bv_copy(&vs->vectors[i]);
        if (!mat[i]) { while (--i >= 0) bv_destroy(mat[i]); free(mat); return 0; }
    }
    int32_t rank = 0, col = 0;
    while (col < d && rank < n) {
        int32_t pivot = -1;
        for (int32_t r = rank; r < n; r++)
            if (bv_get(mat[r], col)) { pivot = r; break; }
        if (pivot < 0) { col++; continue; }
        if (pivot != rank) {
            binary_vector_t *tmp = mat[rank];
            mat[rank] = mat[pivot]; mat[pivot] = tmp;
        }
        for (int32_t r = 0; r < n; r++)
            if (r != rank && bv_get(mat[r], col))
                bv_xor(mat[rank], mat[r], mat[r]);
        rank++; col++;
    }
    for (int32_t i = 0; i < n; i++) bv_destroy(mat[i]);
    free(mat); return rank;
}

typedef struct { int32_t n, d; double density; int64_t orthogonal_pairs; double brute_time, williams_time; } ov_batch_result_t;

int32_t ov_parameter_sweep(const int32_t *n_values, int32_t num_n,
                            double *densities, int32_t num_dens,
                            int32_t d, ov_batch_result_t *results) {
    if (!n_values || !densities || !results) return 0;
    int32_t count = 0;
    for (int32_t ni = 0; ni < num_n; ni++) {
        for (int32_t di = 0; di < num_dens; di++) {
            int32_t n = n_values[ni]; double dens = densities[di];
            ov_instance_t *inst = ov_create_random(n, d, dens, false);
            if (!inst) continue;
            results[count].n = n; results[count].d = d; results[count].density = dens;
            results[count].orthogonal_pairs = ov_count_orthogonal_pairs(inst);
            ov_result_t r = ov_benchmark(inst, ov_brute_force, 1);
            results[count].brute_time = r.time_s;
            r = ov_benchmark(inst, ov_williams, 1);
            results[count].williams_time = r.time_s;
            ov_destroy(inst); count++;
        }
    }
    return count;
}

void ov_print_batch_results(const ov_batch_result_t *results, int32_t count) {
    if (!results || count <= 0) return;
    printf("\n=== OV Parameter Sweep Results ===\n");
    printf("n     | d   | density | orth_pairs | brute_t   | williams_t\n");
    printf("------|-----|---------|------------|-----------|-----------\n");
    for (int32_t i = 0; i < count; i++)
        printf("%-5d | %-3d | %-7.3f | %-10lld | %-9.6f | %-9.6f\n",
               results[i].n, results[i].d, results[i].density,
               (long long)results[i].orthogonal_pairs,
               results[i].brute_time, results[i].williams_time);
}

int32_t ov_min_distance_to_orthogonal(const ov_instance_t *inst) {
    if (!inst || !inst->A || !inst->B) return -1;
    int32_t n = inst->A->num_vectors;
    int32_t min_dist = inst->A->dimension + 1;
    for (int32_t i = 0; i < n; i++)
        for (int32_t j = 0; j < n; j++) {
            int32_t d = bv_dot_product(&inst->A->vectors[i], &inst->B->vectors[j]);
            if (d < min_dist) min_dist = d;
        }
    return min_dist;
}

int32_t ov_max_overlap(const ov_instance_t *inst) {
    if (!inst || !inst->A || !inst->B) return -1;
    int32_t n = inst->A->num_vectors, max_o = 0;
    for (int32_t i = 0; i < n; i++)
        for (int32_t j = 0; j < n; j++) {
            int32_t d = bv_dot_product(&inst->A->vectors[i], &inst->B->vectors[j]);
            if (d > max_o) max_o = d;
        }
    return max_o;
}

void ov_print_implication_chain(void) {
    printf("\n=== Fine-Grained Implication Chain ===\n\n");
    printf("  SETH\n   |\n   | Williams (2005)\n   v\n  OV Conjecture\n   |\n");
    printf("   +--> Pattern Matching [AWY15]\n   +--> Edit Distance [BK15]\n");
    printf("   +--> LCS [ABW15]\n   +--> Frechet [B14]\n   +--> DTW [ABW15]\n");
    printf("   +--> Graph Diameter [RW13]\n   +--> HNN [AW15]\n   +--> Subset Sum [ABHS19]\n\n");
}

double ov_instance_difficulty(const ov_instance_t *inst) {
    if (!inst || !inst->A || !inst->B) return -1.0;
    double score = 0.0;
    double log_n = log2((double)inst->A->num_vectors);
    double ratio = (double)inst->A->dimension / log_n;
    if (ratio > 3.0) score += 4.0;
    else if (ratio > 2.0) score += 3.0;
    else if (ratio > 1.0) score += 2.0;
    else score += 1.0;
    density_stats_t ds = vs_density_stats(inst->A);
    double dev = fabs(ds.mean_density - 0.5);
    if (dev < 0.1) score += 3.0;
    else if (dev < 0.2) score += 2.0;
    else if (dev < 0.3) score += 1.0;
    double ortho_frac = ov_orthogonal_fraction(inst);
    if (ortho_frac < 0.01) score += 3.0;
    else if (ortho_frac < 0.1) score += 2.0;
    else if (ortho_frac < 0.3) score += 1.0;
    return score < 0.0 ? 0.0 : (score > 10.0 ? 10.0 : score);
}

void ov_print_difficulty_report(const ov_instance_t *inst) {
    if (!inst) return;
    double difficulty = ov_instance_difficulty(inst);
    printf("\n=== Instance Difficulty Report ===\n");
    printf("Difficulty score: %.1f / 10.0\n", difficulty);
    printf("Regime: %s\n",
           ov_is_superlogarithmic(inst->A->num_vectors, inst->A->dimension)
           ? "Hard" : "Easier");
    printf("Orthogonal fraction: %.4f\n", ov_orthogonal_fraction(inst));
    density_stats_t ds = vs_density_stats(inst->A);
    printf("Mean density: %.3f\n\n", ds.mean_density);
}

double ov_expected_orthogonal_fraction(int32_t n, int32_t d, double density) {
    double p_overlap = density * density;
    double p_no_overlap = pow(1.0 - p_overlap, (double)d);
    (void)n;
    return p_no_overlap;
}

double ov_expected_orthogonal_count(int32_t n, int32_t d, double density) {
    return ov_expected_orthogonal_fraction(n, d, density) * (double)n * (double)n;
}

void ov_print_expected_vs_actual(const ov_instance_t *inst) {
    if (!inst || !inst->A) return;
    int32_t n = inst->A->num_vectors;
    int32_t d = inst->A->dimension;
    density_stats_t ds = vs_density_stats(inst->A);
    double expected_frac = ov_expected_orthogonal_fraction(n, d, ds.mean_density);
    double actual_frac = ov_orthogonal_fraction(inst);
    int64_t actual_cnt = ov_count_orthogonal_pairs(inst);
    printf("\n=== Expected vs Actual Orthogonal Pairs ===\n");
    printf("Density: %.4f\n", ds.mean_density);
    printf("Expected fraction: %.6f\n", expected_frac);
    printf("Actual fraction:   %.6f\n", actual_frac);
    printf("Actual count:      %lld\n\n", (long long)actual_cnt);
}
