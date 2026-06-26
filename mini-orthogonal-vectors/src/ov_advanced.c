/* ov_advanced.c - Advanced OV Topics (L8)
 * Dynamic OV data structure, Approximate OV,
 * Polynomial method details, Communication complexity simulation.
 */
#include "ov.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ===================================================================
 * L8.4: Dynamic OV Data Structure
 *
 * Maintains sets A and B under insertions and deletions.
 * Supports: query if an orthogonal pair exists.
 *
 * Simple implementation: recompute on each query.
 * Advanced: maintain hash tables for O(1) query at O(nd) update.
 * =================================================================== */

ov_dynamic_t *ov_dynamic_create(int32_t max_n, int32_t d) {
    ov_dynamic_t *dyn = (ov_dynamic_t *)malloc(sizeof(ov_dynamic_t));
    if (!dyn) return NULL;
    dyn->A = vs_create(max_n, d);
    dyn->B = vs_create(max_n, d);
    if (!dyn->A || !dyn->B) {
        vs_destroy(dyn->A); vs_destroy(dyn->B); free(dyn); return NULL;
    }
    dyn->max_a = max_n;
    dyn->max_b = max_n;
    dyn->count_a = 0;
    dyn->count_b = 0;
    dyn->has_answer = false;
    dyn->cached_a = -1;
    dyn->cached_b = -1;
    return dyn;
}

bool ov_dynamic_insert_a(ov_dynamic_t *dyn, const binary_vector_t *v) {
    if (!dyn || !v || dyn->count_a >= dyn->max_a) return false;
    if (v->dimension != dyn->A->dimension) return false;
    bv_copy_into(v, &dyn->A->vectors[dyn->count_a]);
    dyn->A->vectors[dyn->count_a].id = dyn->count_a;
    dyn->count_a++;
    dyn->has_answer = false;  /* invalidate cache */
    return true;
}

bool ov_dynamic_insert_b(ov_dynamic_t *dyn, const binary_vector_t *v) {
    if (!dyn || !v || dyn->count_b >= dyn->max_b) return false;
    if (v->dimension != dyn->B->dimension) return false;
    bv_copy_into(v, &dyn->B->vectors[dyn->count_b]);
    dyn->B->vectors[dyn->count_b].id = dyn->count_b;
    dyn->count_b++;
    dyn->has_answer = false;
    return true;
}

bool ov_dynamic_query(ov_dynamic_t *dyn, int32_t *a_idx, int32_t *b_idx) {
    if (!dyn) return false;
    /* Check cache first */
    if (dyn->has_answer) {
        if (a_idx) *a_idx = dyn->cached_a;
        if (b_idx) *b_idx = dyn->cached_b;
        return true;
    }
    /* Brute force search over active vectors */
    for (int32_t i = 0; i < dyn->count_a; i++) {
        for (int32_t j = 0; j < dyn->count_b; j++) {
            if (bv_are_orthogonal(&dyn->A->vectors[i], &dyn->B->vectors[j])) {
                dyn->has_answer = true;
                dyn->cached_a = i;
                dyn->cached_b = j;
                if (a_idx) *a_idx = i;
                if (b_idx) *b_idx = j;
                return true;
            }
        }
    }
    return false;
}

bool ov_dynamic_delete_a(ov_dynamic_t *dyn, int32_t idx) {
    if (!dyn || idx < 0 || idx >= dyn->count_a) return false;
    /* Lazy deletion: just decrement count; vector at end moves to fill gap */
    if (idx < dyn->count_a - 1)
        bv_copy_into(&dyn->A->vectors[dyn->count_a - 1], &dyn->A->vectors[idx]);
    dyn->count_a--;
    dyn->has_answer = false;
    return true;
}

bool ov_dynamic_delete_b(ov_dynamic_t *dyn, int32_t idx) {
    if (!dyn || idx < 0 || idx >= dyn->count_b) return false;
    if (idx < dyn->count_b - 1)
        bv_copy_into(&dyn->B->vectors[dyn->count_b - 1], &dyn->B->vectors[idx]);
    dyn->count_b--;
    dyn->has_answer = false;
    return true;
}

void ov_dynamic_destroy(ov_dynamic_t *dyn) {
    if (dyn) { vs_destroy(dyn->A); vs_destroy(dyn->B); free(dyn); }
}

/* ===================================================================
 * L8.5: Approximate OV
 *
 * Find a,b such that <a,b> <= epsilon * d (relaxed orthogonality).
 * For epsilon = 0, this is exact OV.
 * For epsilon > 0, this allows up to epsilon*d overlapping 1-bits.
 *
 * This is potentially easier than exact OV for large epsilon.
 * The hardness of approximate OV is an active research question.
 * =================================================================== */

ov_result_t ov_approximate_orthogonal(const ov_instance_t *inst, double epsilon) {
    ov_result_t r = {false, -1, -1, 0, 0.0};
    if (!inst || !inst->A || !inst->B) return r;
    if (epsilon < 0.0) epsilon = 0.0;
    double t0 = ov_get_time_seconds();
    int32_t n = inst->A->num_vectors;
    int32_t d = inst->A->dimension;
    int32_t max_overlap = (int32_t)(epsilon * (double)d);
    if (max_overlap < 0) max_overlap = 0;
    for (int32_t i = 0; i < n; i++) {
        for (int32_t j = 0; j < n; j++) {
            r.ops++;
            int32_t dot = bv_dot_product(&inst->A->vectors[i], &inst->B->vectors[j]);
            if (dot >= 0 && dot <= max_overlap) {
                r.found = true; r.a_index = i; r.b_index = j;
                r.time_s = ov_get_time_seconds() - t0;
                return r;
            }
        }
    }
    r.time_s = ov_get_time_seconds() - t0;
    return r;
}
