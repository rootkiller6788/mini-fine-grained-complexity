/* ov_algorithms.c - OV Algorithms (L5)
 * Brute force, Williams, Meet-in-the-Middle, Four-Russians,
 * Matrix Multiplication, Polynomial Method, Light-Sparse, LSH.
 */
#include "ov.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ---- L5.1: Brute Force O(n^2 * d/32) ---- */

ov_result_t ov_brute_force(const ov_instance_t *inst) {
    ov_result_t r = {false, -1, -1, 0, 0.0};
    if (!inst || !inst->A || !inst->B) return r;
    double t0 = ov_get_time_seconds();
    int32_t n = inst->A->num_vectors;
    for (int32_t i = 0; i < n; i++) {
        for (int32_t j = 0; j < n; j++) {
            r.ops++;
            if (bv_are_orthogonal(&inst->A->vectors[i], &inst->B->vectors[j])) {
                r.found = true; r.a_index = i; r.b_index = j;
                r.time_s = ov_get_time_seconds() - t0;
                return r;
            }
        }
    }
    r.time_s = ov_get_time_seconds() - t0;
    return r;
}

ov_result_t ov_brute_force_early(const ov_instance_t *inst) {
    ov_result_t r = {false, -1, -1, 0, 0.0};
    if (!inst || !inst->A || !inst->B) return r;
    double t0 = ov_get_time_seconds();
    int32_t n = inst->A->num_vectors;
    /* Process in blocks for better cache locality */
    for (int32_t i = 0; i < n; i++) {
        binary_vector_t *va = &inst->A->vectors[i];
        for (int32_t j = 0; j < n; j++) {
            r.ops++;
            /* Early word-level exit: check first word */
            binary_vector_t *vb = &inst->B->vectors[j];
            if ((va->bits[0] & vb->bits[0]) != 0) continue;
            if (bv_are_orthogonal(va, vb)) {
                r.found = true; r.a_index = i; r.b_index = j;
                r.time_s = ov_get_time_seconds() - t0;
                return r;
            }
        }
    }
    r.time_s = ov_get_time_seconds() - t0;
    return r;
}

/* Compare function for sorting by Hamming weight (sparsest first) */
typedef struct { int32_t idx; int32_t weight; } idx_weight_t;
static int cmp_weight(const void *pa, const void *pb) {
    const idx_weight_t *a = (const idx_weight_t *)pa;
    const idx_weight_t *b = (const idx_weight_t *)pb;
    return (a->weight > b->weight) - (a->weight < b->weight);
}

ov_result_t ov_brute_force_sparse_first(const ov_instance_t *inst) {
    ov_result_t r = {false, -1, -1, 0, 0.0};
    if (!inst || !inst->A || !inst->B) return r;
    double t0 = ov_get_time_seconds();
    int32_t n = inst->A->num_vectors;
    /* Sort A indices by Hamming weight (ascending) */
    idx_weight_t *a_order = (idx_weight_t *)malloc((size_t)n * sizeof(idx_weight_t));
    if (!a_order) return ov_brute_force(inst);
    for (int32_t i = 0; i < n; i++) {
        a_order[i].idx = i;
        a_order[i].weight = bv_hamming_weight(&inst->A->vectors[i]);
    }
    qsort(a_order, (size_t)n, sizeof(idx_weight_t), cmp_weight);
    /* Search with sparse-first ordering */
    for (int32_t ii = 0; ii < n; ii++) {
        int32_t i = a_order[ii].idx;
        binary_vector_t *va = &inst->A->vectors[i];
        for (int32_t j = 0; j < n; j++) {
            r.ops++;
            if (bv_are_orthogonal(va, &inst->B->vectors[j])) {
                r.found = true; r.a_index = i; r.b_index = j;
                r.time_s = ov_get_time_seconds() - t0;
                free(a_order); return r;
            }
        }
    }
    r.time_s = ov_get_time_seconds() - t0;
    free(a_order);
    return r;
}

/* ---- L5.2: Williams' Algorithm (Random Projection) ---- */

ov_result_t ov_williams(const ov_instance_t *inst) {
    return ov_williams_configurable(inst, 0, 5);
}

ov_result_t ov_williams_configurable(const ov_instance_t *inst,
                                      int32_t dprime, int32_t num_trials) {
    ov_result_t r = {false, -1, -1, 0, 0.0};
    if (!inst || !inst->A || !inst->B) return r;
    double t0 = ov_get_time_seconds();
    int32_t n = inst->A->num_vectors;
    int32_t d = inst->A->dimension;
    /* Auto-compute dprime if not specified */
    if (dprime <= 0) {
        dprime = (int32_t)(log2((double)n) * 2.0);
        if (dprime < 2) dprime = 2;
        if (dprime > d) dprime = d / 2;
    }
    if (dprime >= d || dprime < 1) return ov_brute_force(inst);
    int32_t *sample = (int32_t *)malloc((size_t)dprime * sizeof(int32_t));
    if (!sample) return ov_brute_force(inst);
    for (int32_t trial = 0; trial < num_trials; trial++) {
        /* Sample dprime random coordinates */
        for (int32_t s = 0; s < dprime; s++)
            sample[s] = rand() % d;
        /* Check all pairs on sampled coordinates first */
        for (int32_t i = 0; i < n; i++) {
            binary_vector_t *va = &inst->A->vectors[i];
            for (int32_t j = 0; j < n; j++) {
                r.ops++;
                binary_vector_t *vb = &inst->B->vectors[j];
                /* Fast check on sampled coordinates */
                bool possible = true;
                for (int32_t s = 0; s < dprime && possible; s++)
                    if (bv_get(va, sample[s]) && bv_get(vb, sample[s]))
                        possible = false;
                if (possible && bv_are_orthogonal(va, vb)) {
                    r.found = true; r.a_index = i; r.b_index = j;
                    r.time_s = ov_get_time_seconds() - t0;
                    free(sample); return r;
                }
            }
        }
    }
    r.time_s = ov_get_time_seconds() - t0;
    free(sample);
    return r;
}

/* ---- L5.3: Meet-in-the-Middle O(n * 2^{d/2}) ---- */

ov_result_t ov_meet_in_middle(const ov_instance_t *inst) {
    ov_result_t r = {false, -1, -1, 0, 0.0};
    if (!inst || !inst->A || !inst->B) return r;
    int32_t d = inst->A->dimension;
    /* Only practical for small d (<= 24 to keep 2^{d/2} <= 4096) */
    if (d > 24) return ov_brute_force(inst);
    double t0 = ov_get_time_seconds();
    int32_t n = inst->A->num_vectors;
    int32_t half = d / 2;
    int32_t half2 = d - half;
    int32_t size = 1 << half;
    /* For each possible left half pattern, list B vectors matching it */
    /* Use a simple array of linked lists */
    typedef struct node_t { int32_t idx; struct node_t *next; } node_t;
    node_t **buckets = (node_t **)calloc((size_t)size, sizeof(node_t *));
    if (!buckets) return ov_brute_force(inst);
    /* Bucket B vectors by their left-half bit pattern */
    for (int32_t j = 0; j < n; j++) {
        r.ops++;
        int32_t pattern = 0;
        for (int32_t b = 0; b < half; b++)
            if (bv_get(&inst->B->vectors[j], b))
                pattern |= (1 << b);
        node_t *nd = (node_t *)malloc(sizeof(node_t));
        if (!nd) continue;
        nd->idx = j; nd->next = buckets[pattern];
        buckets[pattern] = nd;
    }
    /* For each A vector, check B vectors with complementary left half */
    for (int32_t i = 0; i < n; i++) {
        int32_t a_pattern = 0;
        for (int32_t b = 0; b < half; b++)
            if (bv_get(&inst->A->vectors[i], b))
                a_pattern |= (1 << b);
        /* B vectors must have 0s where A has 1s in left half */
        /* Any B pattern that is subset of ~a_pattern is candidate */
        int32_t complement = ((1 << half) - 1) & ~a_pattern;
        /* Check all subsets of complement (using submask enumeration) */
        int32_t sub = complement;
        do {
            for (node_t *nd = buckets[sub]; nd; nd = nd->next) {
                r.ops++;
                if (bv_are_orthogonal(&inst->A->vectors[i],
                                       &inst->B->vectors[nd->idx])) {
                    r.found = true; r.a_index = i; r.b_index = nd->idx;
                    r.time_s = ov_get_time_seconds() - t0;
                    /* Cleanup */
                    for (int32_t s = 0; s < size; s++) {
                        node_t *cur = buckets[s];
                        while (cur) { node_t *nx = cur->next; free(cur); cur = nx; }
                    }
                    free(buckets); return r;
                }
            }
            sub = (sub - 1) & complement;
        } while (sub != complement);
    }
    /* Cleanup */
    for (int32_t s = 0; s < size; s++) {
        node_t *cur = buckets[s];
        while (cur) { node_t *nx = cur->next; free(cur); cur = nx; }
    }
    free(buckets);
    r.time_s = ov_get_time_seconds() - t0;
    return r;
}

/* ---- L5.4: Four-Russians O(n^2 * d / log n) ---- */

ov_result_t ov_four_russians(const ov_instance_t *inst) {
    ov_result_t r = {false, -1, -1, 0, 0.0};
    if (!inst || !inst->A || !inst->B) return r;
    double t0 = ov_get_time_seconds();
    int32_t n = inst->A->num_vectors;
    int32_t d = inst->A->dimension;
    /* Block size: b = floor(log2(n) / 2) */
    int32_t b = (int32_t)(log2((double)n) / 2.0);
    if (b < 2) b = 2;
    if (b > 16) b = 16; /* 2^16 = 65536 precomputed entries, manageable */
    int32_t num_blocks = (d + b - 1) / b;
    int32_t table_size = 1 << b;
    /* Precompute orthogonality for all 2^b x 2^b block patterns */
    /* orthog_table[x][y] = true if x AND y == 0 */
    bool *orthog_table = (bool *)malloc((size_t)table_size * (size_t)table_size * sizeof(bool));
    if (!orthog_table) return ov_brute_force(inst);
    for (int32_t x = 0; x < table_size; x++)
        for (int32_t y = 0; y < table_size; y++)
            orthog_table[x * table_size + y] = ((x & y) == 0);
    /* Extract block values for all vectors */
    int32_t **A_blocks = (int32_t **)malloc((size_t)n * sizeof(int32_t *));
    int32_t **B_blocks = (int32_t **)malloc((size_t)n * sizeof(int32_t *));
    if (!A_blocks || !B_blocks) {
        free(orthog_table); free(A_blocks); free(B_blocks);
        return ov_brute_force(inst);
    }
    for (int32_t i = 0; i < n; i++) {
        A_blocks[i] = (int32_t *)malloc((size_t)num_blocks * sizeof(int32_t));
        B_blocks[i] = (int32_t *)malloc((size_t)num_blocks * sizeof(int32_t));
        if (!A_blocks[i] || !B_blocks[i]) {
            /* Partial cleanup on failure */
            for (int32_t j = 0; j <= i; j++) { free(A_blocks[j]); free(B_blocks[j]); }
            free(A_blocks); free(B_blocks); free(orthog_table);
            return ov_brute_force(inst);
        }
    }
    /* Extract blocks */
    for (int32_t i = 0; i < n; i++) {
        for (int32_t blk = 0; blk < num_blocks; blk++) {
            int32_t aval = 0, bval = 0;
            for (int32_t bit = 0; bit < b; bit++) {
                int32_t pos = blk * b + bit;
                if (pos < d) {
                    if (bv_get(&inst->A->vectors[i], pos)) aval |= (1 << bit);
                    if (bv_get(&inst->B->vectors[i], pos)) bval |= (1 << bit);
                }
            }
            A_blocks[i][blk] = aval;
            B_blocks[i][blk] = bval;
        }
    }
    /* Four-Russians lookup: check each pair using precomputed table */
    for (int32_t i = 0; i < n; i++) {
        for (int32_t j = 0; j < n; j++) {
            r.ops++;
            bool orthogonal = true;
            for (int32_t blk = 0; blk < num_blocks && orthogonal; blk++) {
                if (!orthog_table[A_blocks[i][blk] * table_size + B_blocks[j][blk]])
                    orthogonal = false;
            }
            if (orthogonal) {
                r.found = true; r.a_index = i; r.b_index = j;
                r.time_s = ov_get_time_seconds() - t0;
                goto cleanup;
            }
        }
    }
cleanup:
    r.time_s = ov_get_time_seconds() - t0;
    for (int32_t i = 0; i < n; i++) { free(A_blocks[i]); free(B_blocks[i]); }
    free(A_blocks); free(B_blocks); free(orthog_table);
    return r;
}

/* ---- L5.5: Matrix Multiplication Approach ---- */

ov_result_t ov_via_matrix_mult(const ov_instance_t *inst) {
    ov_result_t r = {false, -1, -1, 0, 0.0};
    if (!inst || !inst->A || !inst->B) return r;
    double t0 = ov_get_time_seconds();
    int32_t n = inst->A->num_vectors;
    int32_t d = inst->A->dimension;
    /* Compute n x n result matrix P = A * B^T
     * P[i][j] = sum_k A[i][k] * B[j][k] = <A[i], B[j]>
     * Using standard O(n^2*d) multiply (educational; fast MM is complex) */
    int32_t *P = (int32_t *)calloc((size_t)n * (size_t)n, sizeof(int32_t));
    if (!P) return ov_brute_force(inst);
    for (int32_t i = 0; i < n; i++) {
        for (int32_t j = 0; j < n; j++) {
            r.ops++;
            int32_t dot = bv_dot_product(&inst->A->vectors[i], &inst->B->vectors[j]);
            P[i * n + j] = dot;
            if (dot == 0) {
                r.found = true; r.a_index = i; r.b_index = j;
                r.time_s = ov_get_time_seconds() - t0;
                free(P); return r;
            }
        }
    }
    r.time_s = ov_get_time_seconds() - t0;
    free(P);
    return r;
}

/* ---- L5.6: Polynomial Method ---- */

ov_result_t ov_polynomial_method(const ov_instance_t *inst) {
    /* The polynomial method evaluates P(a,b) = prod(1 - a_i*b_i)
     * over GF(2). For small d, use meet-in-the-middle on polynomial.
     * This implementation is a simplified version for educational purposes. */
    ov_result_t r = {false, -1, -1, 0, 0.0};
    if (!inst || !inst->A || !inst->B) return r;
    /* For small d (<= 20), this reduces to subset-sum-like approach */
    int32_t d = inst->A->dimension;
    if (d <= 20) return ov_meet_in_middle(inst);
    /* For larger d, the polynomial method requires sophisticated
     * fast multipoint evaluation which we simplify to brute force. */
    return ov_brute_force(inst);
}

/* ---- L5.7: Light-Sparse Vectors ---- */

ov_result_t ov_light_sparse(const ov_instance_t *inst) {
    ov_result_t r = {false, -1, -1, 0, 0.0};
    if (!inst || !inst->A || !inst->B) return r;
    /* Check if vectors are sufficiently sparse for sparse representation */
    density_stats_t dsA = vs_density_stats(inst->A);
    density_stats_t dsB = vs_density_stats(inst->B);
    double avg_density = (dsA.mean_density + dsB.mean_density) / 2.0;
    if (avg_density > 0.1) return ov_brute_force(inst);
    double t0 = ov_get_time_seconds();
    int32_t n = inst->A->num_vectors;
    int32_t d = inst->A->dimension;
    /* Build sparse index: for each vector, list positions of 1-bits */
    int32_t max_ones = (int32_t)((double)d * 0.3 + 1);
    int32_t **sparse_A = (int32_t **)malloc((size_t)n * sizeof(int32_t *));
    int32_t **sparse_B = (int32_t **)malloc((size_t)n * sizeof(int32_t *));
    int32_t *count_A = (int32_t *)calloc((size_t)n, sizeof(int32_t));
    int32_t *count_B = (int32_t *)calloc((size_t)n, sizeof(int32_t));
    if (!sparse_A || !sparse_B || !count_A || !count_B) {
        free(sparse_A); free(sparse_B); free(count_A); free(count_B);
        return ov_brute_force(inst);
    }
    for (int32_t i = 0; i < n; i++) {
        sparse_A[i] = (int32_t *)malloc((size_t)max_ones * sizeof(int32_t));
        sparse_B[i] = (int32_t *)malloc((size_t)max_ones * sizeof(int32_t));
        if (!sparse_A[i] || !sparse_B[i]) { r = ov_brute_force(inst); goto spclean; }
        for (int32_t k = 0; k < d; k++) {
            if (bv_get(&inst->A->vectors[i], k) && count_A[i] < max_ones)
                sparse_A[i][count_A[i]++] = k;
            if (bv_get(&inst->B->vectors[i], k) && count_B[i] < max_ones)
                sparse_B[i][count_B[i]++] = k;
        }
    }
    /* Check orthogonality using sparse lists: intersection of one-positions */
    for (int32_t i = 0; i < n; i++) {
        for (int32_t j = 0; j < n; j++) {
            r.ops++;
            /* Check if any 1-bit position of A[i] also appears in B[j] */
            bool orth = true;
            for (int32_t ka = 0; ka < count_A[i] && orth; ka++) {
                int32_t pos = sparse_A[i][ka];
                for (int32_t kb = 0; kb < count_B[j]; kb++)
                    if (sparse_B[j][kb] == pos) { orth = false; break; }
            }
            if (orth) {
                r.found = true; r.a_index = i; r.b_index = j;
                r.time_s = ov_get_time_seconds() - t0;
                goto spclean;
            }
        }
    }
spclean:
    r.time_s = ov_get_time_seconds() - t0;
    for (int32_t i = 0; i < n; i++) { free(sparse_A[i]); free(sparse_B[i]); }
    free(sparse_A); free(sparse_B); free(count_A); free(count_B);
    return r;
}

/* ---- L5.8: Locality-Sensitive Hashing ---- */

ov_result_t ov_lsh(const ov_instance_t *inst, int32_t L, int32_t k) {
    ov_result_t r = {false, -1, -1, 0, 0.0};
    if (!inst || !inst->A || !inst->B) return r;
    if (L <= 0) L = 5;
    if (k <= 0) k = 1;
    double t0 = ov_get_time_seconds();
    int32_t n = inst->A->num_vectors;
    int32_t d = inst->A->dimension;
    /* For each hash table, sample k random bit positions as hash key */
    for (int32_t table = 0; table < L; table++) {
        int32_t *hash_bits = (int32_t *)malloc((size_t)k * sizeof(int32_t));
        if (!hash_bits) continue;
        for (int32_t h = 0; h < k; h++)
            hash_bits[h] = rand() % d;
        /* Hash all A and B vectors; bucket by hash value */
        for (int32_t i = 0; i < n; i++) {
            for (int32_t j = 0; j < n; j++) {
                r.ops++;
                bool all_zero = true;
                for (int32_t h = 0; h < k && all_zero; h++) {
                    int32_t p = hash_bits[h];
                    if (bv_get(&inst->A->vectors[i], p) &&
                        bv_get(&inst->B->vectors[j], p))
                        all_zero = false;
                }
                if (all_zero && bv_are_orthogonal(&inst->A->vectors[i],
                                                    &inst->B->vectors[j])) {
                    r.found = true; r.a_index = i; r.b_index = j;
                    r.time_s = ov_get_time_seconds() - t0;
                    free(hash_bits); return r;
                }
            }
        }
        free(hash_bits);
    }
    r.time_s = ov_get_time_seconds() - t0;
    return r;
}

/* ---- L5.9: Counting Variants ---- */

int64_t ov_count_orthogonal_pairs(const ov_instance_t *inst) {
    if (!inst || !inst->A || !inst->B) return 0;
    int32_t n = inst->A->num_vectors;
    int64_t count = 0;
    for (int32_t i = 0; i < n; i++)
        for (int32_t j = 0; j < n; j++)
            if (bv_are_orthogonal(&inst->A->vectors[i], &inst->B->vectors[j]))
                count++;
    return count;
}

bool ov_has_k_orthogonal_pairs(const ov_instance_t *inst, int64_t K) {
    if (!inst || !inst->A || !inst->B) return (K <= 0);
    int32_t n = inst->A->num_vectors;
    int64_t count = 0;
    for (int32_t i = 0; i < n; i++) {
        for (int32_t j = 0; j < n; j++) {
            if (bv_are_orthogonal(&inst->A->vectors[i], &inst->B->vectors[j])) {
                count++;
                if (count >= K) return true;
            }
        }
    }
    return count >= K;
}

int64_t ov_find_all_orthogonal_pairs(const ov_instance_t *inst,
                                      int32_t *ai, int32_t *bi, int64_t max) {
    if (!inst || !inst->A || !inst->B || !ai || !bi) return 0;
    int32_t n = inst->A->num_vectors;
    int64_t found = 0;
    for (int32_t i = 0; i < n && found < max; i++) {
        for (int32_t j = 0; j < n && found < max; j++) {
            if (bv_are_orthogonal(&inst->A->vectors[i], &inst->B->vectors[j])) {
                ai[found] = i;
                bi[found] = j;
                found++;
            }
        }
    }
    return found;
}

/* ---- Streaming Online Algorithm ---- */

ov_result_t ov_streaming_online(const ov_instance_t *inst) {
    /* Streaming model: A vectors arrive online.
     * For each arriving a, check against all stored B vectors.
     * This is just brute force but with the "streaming" abstraction. */
    return ov_brute_force(inst);
}
