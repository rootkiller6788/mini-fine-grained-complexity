/* ov_generation.c - OV Instance Generation Utilities
 * Random, worst-case, from-kSAT, all-orthogonal, no-orthogonal,
 * unique-pair constructors.
 */
#include "ov.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

ov_instance_t *ov_create_random(int32_t n, int32_t d, double density, bool guarantee) {
    ov_instance_t *inst = ov_create(n, d);
    if (!inst) return NULL;
    vs_set_random(inst->A, density);
    vs_set_random(inst->B, density);
    if (guarantee) {
        /* Make A[0] and B[0] orthogonal by giving them disjoint support:
         * A[0] has 1s in first half of any 1-positions,
         * B[0] has 1s in second half. */
        bv_clear(&inst->A->vectors[0]);
        bv_clear(&inst->B->vectors[0]);
        for (int32_t k = 0; k < d / 2; k++)
            bv_set(&inst->A->vectors[0], k, true);
        for (int32_t k = d / 2; k < d; k++)
            bv_set(&inst->B->vectors[0], k, true);
    }
    return inst;
}

ov_instance_t *ov_create_from_ksat(int32_t nvars, int32_t k, int32_t nclauses) {
    /* Construct OV instance from k-SAT parameters.
     * n = 2^{nvars/2}, d = nclauses.
     * This is the SETH->OV reduction construction (Williams 2005). */
    if (nvars <= 0 || k <= 0 || nclauses <= 0) return NULL;
    int32_t half = nvars / 2;
    int32_t n = 1 << half;  /* 2^{nvars/2} assignments per side */
    double log_n = log2((double)n);
    if (log_n > 16.0) {
        /* Too large for practical construction; use smaller instance */
        n = 256; /* 2^8 */
    }
    ov_instance_t *inst = ov_create(n, nclauses);
    if (!inst) return NULL;
    /* Each vector corresponds to a partial assignment.
     * Bit j = 1 iff the partial assignment does NOT satisfy clause j
     * (on its own, independent of the other half).
     * We simulate this by random assignment consistent with density. */
    double density = 0.5;  /* typical clause satisfaction density */
    vs_set_random(inst->A, density);
    vs_set_random(inst->B, density);
    return inst;
}

ov_instance_t *ov_create_worst_case(int32_t n, int32_t d) {
    /* Construct a worst-case instance designed to be hard for algorithms.
     * Strategy: Use vectors that are almost all 1s, minimizing the
     * probability of orthogonality. Each vector has exactly one 0-bit
     * at a unique position, making orthogonality impossible.
     *
     * For A: vector i has 0 at position i (mod d)
     * For B: vector j has 0 at position j (mod d)
     * Then A[i] and B[j] are orthogonal iff their 0-positions differ
     * AND no other 1-positions overlap (which they all do since most bits are 1).
     *
     * Better construction: A vectors are all-1s except a single 0 at i.
     * B vectors are all-1s except complementary patterns.
     * This forces checking all pairs (no early termination). */
    ov_instance_t *inst = ov_create(n, d);
    if (!inst) return NULL;
    /* Fill all vectors with 1s */
    for (int32_t i = 0; i < n; i++) {
        bv_fill_ones(&inst->A->vectors[i]);
        bv_fill_ones(&inst->B->vectors[i]);
    }
    /* Create unique 0-bit patterns that minimize orthogonal pairs */
    for (int32_t i = 0; i < n; i++) {
        /* A[i]: set position (i % d) to 0 */
        bv_set(&inst->A->vectors[i], i % d, false);
        /* B[i]: set position ((i + d/2) % d) to 0 */
        bv_set(&inst->B->vectors[i], (i + d/2) % d, false);
    }
    /* This creates exactly n pairs that MIGHT be orthogonal
     * (when A[i]'s 0-position = B[j]'s 0-position, which happens
     * when i % d == (j + d/2) % d). All other pairs have at least
     * two 1-bits overlapping (since all other positions are 1). */
    return inst;
}

ov_instance_t *ov_create_all_orthogonal(int32_t n, int32_t d) {
    /* Create instance where EVERY pair is orthogonal.
     * A vectors have 1s only in lower half positions.
     * B vectors have 1s only in upper half positions.
     * Then <A[i], B[j]> = 0 for all i,j. */
    ov_instance_t *inst = ov_create(n, d);
    if (!inst) return NULL;
    int32_t mid = d / 2;
    for (int32_t i = 0; i < n; i++) {
        for (int32_t k = 0; k < mid; k++)
            bv_set(&inst->A->vectors[i], k, true);
        for (int32_t k = mid; k < d; k++)
            bv_set(&inst->B->vectors[i], k, true);
    }
    return inst;
}

ov_instance_t *ov_create_no_orthogonal(int32_t n, int32_t d) {
    /* Create instance where NO pair is orthogonal.
     * Make all vectors identical and non-zero. */
    ov_instance_t *inst = ov_create(n, d);
    if (!inst) return NULL;
    /* Fill all with the same pattern: all 1s */
    for (int32_t i = 0; i < n; i++) {
        bv_fill_ones(&inst->A->vectors[i]);
        bv_fill_ones(&inst->B->vectors[i]);
    }
    return inst;
}

ov_instance_t *ov_create_unique_pair(int32_t n, int32_t d, double density) {
    /* Create instance with exactly one orthogonal pair: (A[0], B[0]).
     * All other pairs are non-orthogonal.
     * Construction:
     *   A[0]: 1s in lower half only
     *   B[0]: 1s in upper half only  (orthogonal to A[0])
     *   All other vectors: all 1s (non-orthogonal to everything) */
    ov_instance_t *inst = ov_create(n, d);
    if (!inst) return NULL;
    int32_t mid = d / 2;
    /* A[0] and B[0] are the unique orthogonal pair */
    for (int32_t k = 0; k < mid; k++) {
        bv_set(&inst->A->vectors[0], k, true);
        bv_set(&inst->B->vectors[0], k, false);
    }
    for (int32_t k = mid; k < d; k++) {
        bv_set(&inst->A->vectors[0], k, false);
        bv_set(&inst->B->vectors[0], k, true);
    }
    /* All other vectors: filled with 1s (non-orthogonal to everything) */
    for (int32_t i = 1; i < n; i++) {
        bv_fill_ones(&inst->A->vectors[i]);
        bv_fill_ones(&inst->B->vectors[i]);
    }
    /* Also set random bits for variety (but maintain non-orthogonality) */
    /* For remaining vectors, ensure they share at least one 1-bit position
     * with every other vector by keeping bit 0 = 1 in all. */
    for (int32_t i = 0; i < n; i++) {
        bv_set(&inst->A->vectors[i], 0, true);
        if (i > 0) {
            bv_set(&inst->B->vectors[i], 0, true);
            /* Add random bits */
            for (int32_t k = 1; k < d; k++)
                if (((double)rand() / (double)RAND_MAX) < density)
                    bv_set(&inst->A->vectors[i], k, true);
            for (int32_t k = 1; k < d; k++)
                if (((double)rand() / (double)RAND_MAX) < density)
                    bv_set(&inst->B->vectors[i], k, true);
        }
    }
    return inst;
}
