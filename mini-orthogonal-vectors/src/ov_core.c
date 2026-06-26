/* ov_core.c - Core OV data structures and operations (L1-L3)
 * Implements binary vectors, vector sets, OV instances, packed bit ops.
 */
#include "ov.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>

#ifdef _MSC_VER
#include <intrin.h>
#define BV_POPCOUNT32 __popcnt
#else
#define BV_POPCOUNT32 __builtin_popcount
#endif

/* ===================================================================
 * L1: Binary Vector Creation, Destruction, Copying
 * =================================================================== */

binary_vector_t *bv_create(int32_t dim) {
    if (dim <= 0) return NULL;
    binary_vector_t *v = (binary_vector_t *)malloc(sizeof(binary_vector_t));
    if (!v) return NULL;
    v->dimension = dim;
    v->id = 0;
    int32_t words = (dim + 31) / 32;
    v->bits = (uint32_t *)calloc((size_t)words, sizeof(uint32_t));
    if (!v->bits) { free(v); return NULL; }
    return v;
}

binary_vector_t *bv_copy(const binary_vector_t *src) {
    if (!src) return NULL;
    binary_vector_t *v = bv_create(src->dimension);
    if (!v) return NULL;
    int32_t words = (src->dimension + 31) / 32;
    memcpy(v->bits, src->bits, (size_t)words * sizeof(uint32_t));
    v->id = src->id;
    return v;
}

void bv_copy_into(const binary_vector_t *src, binary_vector_t *dst) {
    if (!src || !dst || src->dimension != dst->dimension) return;
    int32_t words = (src->dimension + 31) / 32;
    memcpy(dst->bits, src->bits, (size_t)words * sizeof(uint32_t));
    dst->id = src->id;
}

void bv_destroy(binary_vector_t *v) {
    if (v) { free(v->bits); free(v); }
}

/* ===================================================================
 * L1: Bit Access and Manipulation
 * =================================================================== */

bool bv_get(const binary_vector_t *v, int32_t pos) {
    if (!v || pos < 0 || pos >= v->dimension) return false;
    return (v->bits[(uint32_t)pos / 32] >> ((uint32_t)pos % 32)) & 1U;
}

void bv_set(binary_vector_t *v, int32_t pos, bool val) {
    if (!v || pos < 0 || pos >= v->dimension) return;
    uint32_t idx = (uint32_t)pos / 32;
    uint32_t mask = 1U << ((uint32_t)pos % 32);
    if (val) v->bits[idx] |= mask;
    else     v->bits[idx] &= ~mask;
}

void bv_flip(binary_vector_t *v, int32_t pos) {
    if (!v || pos < 0 || pos >= v->dimension) return;
    v->bits[(uint32_t)pos / 32] ^= (1U << ((uint32_t)pos % 32));
}

void bv_clear(binary_vector_t *v) {
    if (!v) return;
    int32_t words = (v->dimension + 31) / 32;
    memset(v->bits, 0, (size_t)words * sizeof(uint32_t));
}

void bv_fill_ones(binary_vector_t *v) {
    if (!v) return;
    int32_t words = (v->dimension + 31) / 32;
    memset(v->bits, 0xFF, (size_t)words * sizeof(uint32_t));
    /* Mask the last word to avoid spurious bits beyond dimension */
    int32_t rem = v->dimension % 32;
    if (rem > 0) v->bits[words - 1] &= (1U << rem) - 1;
}

void bv_set_random(binary_vector_t *v, double density) {
    if (!v) return;
    for (int32_t i = 0; i < v->dimension; i++)
        bv_set(v, i, ((double)rand() / (double)RAND_MAX) < density);
}

bool bv_set_from_string(binary_vector_t *v, const char *s) {
    if (!v || !s) return false;
    int32_t len = (int32_t)strlen(s);
    if (len != v->dimension) return false;
    for (int32_t i = 0; i < len; i++) {
        if (s[i] == '1') bv_set(v, i, true);
        else if (s[i] == '0') bv_set(v, i, false);
        else return false;
    }
    return true;
}

/* ===================================================================
 * L1: Dot Product and Orthogonality
 * =================================================================== */

int32_t bv_dot_product(const binary_vector_t *a, const binary_vector_t *b) {
    if (!a || !b || a->dimension != b->dimension) return -1;
    int32_t words = (a->dimension + 31) / 32;
    int32_t sum = 0;
    for (int32_t i = 0; i < words; i++) {
        uint32_t w = a->bits[i] & b->bits[i];
        sum += BV_POPCOUNT32(w);
    }
    return sum;
}

bool bv_are_orthogonal(const binary_vector_t *a, const binary_vector_t *b) {
    if (!a || !b || a->dimension != b->dimension) return false;
    int32_t words = (a->dimension + 31) / 32;
    for (int32_t i = 0; i < words; i++) {
        if ((a->bits[i] & b->bits[i]) != 0) return false;
    }
    return true;
}

int32_t bv_hamming_weight(const binary_vector_t *v) {
    if (!v) return 0;
    int32_t words = (v->dimension + 31) / 32;
    int32_t w = 0;
    for (int32_t i = 0; i < words; i++)
        w += BV_POPCOUNT32(v->bits[i]);
    return w;
}

int32_t bv_hamming_distance(const binary_vector_t *a, const binary_vector_t *b) {
    if (!a || !b || a->dimension != b->dimension) return -1;
    int32_t words = (a->dimension + 31) / 32;
    int32_t dist = 0;
    for (int32_t i = 0; i < words; i++)
        dist += BV_POPCOUNT32(a->bits[i] ^ b->bits[i]);
    return dist;
}

bool bv_orthogonal_on_subset(const binary_vector_t *a, const binary_vector_t *b,
                              const int32_t *pos, int32_t n) {
    if (!a || !b || !pos) return false;
    for (int32_t k = 0; k < n; k++) {
        int32_t p = pos[k];
        if (p < 0 || p >= a->dimension || p >= b->dimension) continue;
        if (bv_get(a, p) && bv_get(b, p)) return false;
    }
    return true;
}

/* ===================================================================
 * L1: Similarity Metrics
 * =================================================================== */

double bv_jaccard(const binary_vector_t *a, const binary_vector_t *b) {
    if (!a || !b || a->dimension != b->dimension) return -1.0;
    int32_t words = (a->dimension + 31) / 32;
    int32_t intersection = 0, union_count = 0;
    for (int32_t i = 0; i < words; i++) {
        uint32_t aw = a->bits[i], bw = b->bits[i];
        intersection += BV_POPCOUNT32(aw & bw);
        union_count  += BV_POPCOUNT32(aw | bw);
    }
    if (union_count == 0) return 1.0;
    return (double)intersection / (double)union_count;
}

double bv_cosine(const binary_vector_t *a, const binary_vector_t *b) {
    if (!a || !b) return -1.0;
    int32_t dot = bv_dot_product(a, b);
    if (dot < 0) return -1.0;
    int32_t wa = bv_hamming_weight(a);
    int32_t wb = bv_hamming_weight(b);
    if (wa == 0 || wb == 0) return 0.0;
    return (double)dot / sqrt((double)wa * (double)wb);
}

/* ===================================================================
 * L1: I/O
 * =================================================================== */

void bv_print(const binary_vector_t *v) {
    bv_fprint(stdout, v);
}

void bv_fprint(FILE *fp, const binary_vector_t *v) {
    if (!fp || !v) return;
    fprintf(fp, "V[%d] d=%d: [", v->id, v->dimension);
    int32_t lim = v->dimension < 64 ? v->dimension : 64;
    for (int32_t i = 0; i < lim; i++)
        fputc(bv_get(v, i) ? '1' : '0', fp);
    if (v->dimension > 64) fprintf(fp, "...(%d more)", v->dimension - 64);
    fprintf(fp, "]\n");
}

/* ===================================================================
 * L1: Vector Set Operations
 * =================================================================== */

vector_set_t *vs_create(int32_t n, int32_t dim) {
    if (n <= 0 || dim <= 0) return NULL;
    vector_set_t *vs = (vector_set_t *)malloc(sizeof(vector_set_t));
    if (!vs) return NULL;
    vs->num_vectors = n;
    vs->dimension = dim;
    vs->vectors = (binary_vector_t *)calloc((size_t)n, sizeof(binary_vector_t));
    if (!vs->vectors) { free(vs); return NULL; }
    for (int32_t i = 0; i < n; i++) {
        vs->vectors[i].dimension = dim;
        vs->vectors[i].id = i;
        int32_t words = (dim + 31) / 32;
        vs->vectors[i].bits = (uint32_t *)calloc((size_t)words, sizeof(uint32_t));
        if (!vs->vectors[i].bits) {
            for (int32_t j = 0; j < i; j++) free(vs->vectors[j].bits);
            free(vs->vectors); free(vs);
            return NULL;
        }
    }
    return vs;
}

vector_set_t *vs_copy(const vector_set_t *src) {
    if (!src) return NULL;
    vector_set_t *vs = vs_create(src->num_vectors, src->dimension);
    if (!vs) return NULL;
    int32_t words = (src->dimension + 31) / 32;
    for (int32_t i = 0; i < src->num_vectors; i++)
        memcpy(vs->vectors[i].bits, src->vectors[i].bits,
               (size_t)words * sizeof(uint32_t));
    return vs;
}

void vs_destroy(vector_set_t *vs) {
    if (!vs) return;
    for (int32_t i = 0; i < vs->num_vectors; i++)
        free(vs->vectors[i].bits);
    free(vs->vectors);
    free(vs);
}

binary_vector_t *vs_get(vector_set_t *vs, int32_t i) {
    if (!vs || i < 0 || i >= vs->num_vectors) return NULL;
    return &vs->vectors[i];
}

void vs_set_random(vector_set_t *vs, double density) {
    if (!vs) return;
    for (int32_t i = 0; i < vs->num_vectors; i++)
        bv_set_random(&vs->vectors[i], density);
}

int32_t vs_set_from_strings(vector_set_t *vs, const char **pats, int32_t n) {
    if (!vs || !pats || n != vs->num_vectors) return 0;
    int32_t ok = 0;
    for (int32_t i = 0; i < n; i++)
        if (bv_set_from_string(&vs->vectors[i], pats[i])) ok++;
    return ok;
}

density_stats_t vs_density_stats(const vector_set_t *vs) {
    density_stats_t ds;
    memset(&ds, 0, sizeof(ds));
    if (!vs || vs->num_vectors == 0) return ds;
    ds.num_vectors = vs->num_vectors;
    double total_density = 0.0;
    ds.min_density = 1.0;
    ds.max_density = 0.0;
    for (int32_t i = 0; i < vs->num_vectors; i++) {
        int32_t w = bv_hamming_weight(&vs->vectors[i]);
        double dens = (double)w / (double)vs->dimension;
        total_density += dens;
        if (dens < ds.min_density) ds.min_density = dens;
        if (dens > ds.max_density) ds.max_density = dens;
        if (w == 0) ds.num_zero_vectors++;
        if (w == vs->dimension) ds.num_ones_vectors++;
    }
    ds.mean_density = total_density / (double)vs->num_vectors;
    /* Compute stddev */
    double sum_sq = 0.0;
    for (int32_t i = 0; i < vs->num_vectors; i++) {
        double w = (double)bv_hamming_weight(&vs->vectors[i]);
        double dens = w / (double)vs->dimension;
        double diff = dens - ds.mean_density;
        sum_sq += diff * diff;
    }
    ds.stddev_density = sqrt(sum_sq / (double)vs->num_vectors);
    return ds;
}

void vs_print(const vector_set_t *vs) {
    if (!vs) return;
    printf("Vector Set: %d vectors, dim=%d\n", vs->num_vectors, vs->dimension);
    for (int32_t i = 0; i < vs->num_vectors && i < 16; i++)
        bv_print(&vs->vectors[i]);
    if (vs->num_vectors > 16)
        printf("  ... (%d more vectors)\n", vs->num_vectors - 16);
}

/* ===================================================================
 * L1: OV Instance Operations
 * =================================================================== */

ov_instance_t *ov_create(int32_t n, int32_t d) {
    ov_instance_t *inst = (ov_instance_t *)malloc(sizeof(ov_instance_t));
    if (!inst) return NULL;
    inst->A = vs_create(n, d);
    inst->B = vs_create(n, d);
    if (!inst->A || !inst->B) {
        vs_destroy(inst->A); vs_destroy(inst->B); free(inst);
        return NULL;
    }
    return inst;
}

ov_instance_t *ov_copy(const ov_instance_t *inst) {
    if (!inst) return NULL;
    ov_instance_t *cpy = (ov_instance_t *)malloc(sizeof(ov_instance_t));
    if (!cpy) return NULL;
    cpy->A = vs_copy(inst->A);
    cpy->B = vs_copy(inst->B);
    if (!cpy->A || !cpy->B) {
        vs_destroy(cpy->A); vs_destroy(cpy->B); free(cpy);
        return NULL;
    }
    return cpy;
}

void ov_destroy(ov_instance_t *inst) {
    if (!inst) return;
    vs_destroy(inst->A);
    vs_destroy(inst->B);
    free(inst);
}

int64_t ov_num_pairs(const ov_instance_t *inst) {
    if (!inst || !inst->A || !inst->B) return 0;
    return (int64_t)inst->A->num_vectors * (int64_t)inst->B->num_vectors;
}

/* ===================================================================
 * L3: Packed Bit Operations
 * =================================================================== */

int32_t bv_num_words(int32_t dim) {
    return (dim + 31) / 32;
}

void bv_word_bit_index(int32_t pos, int32_t *wi, uint32_t *bm) {
    *wi = pos / 32;
    *bm = 1U << ((uint32_t)pos % 32);
}

void bv_and(const binary_vector_t *a, const binary_vector_t *b, binary_vector_t *r) {
    if (!a || !b || !r) return;
    if (a->dimension != b->dimension || a->dimension != r->dimension) return;
    int32_t words = (a->dimension + 31) / 32;
    for (int32_t i = 0; i < words; i++)
        r->bits[i] = a->bits[i] & b->bits[i];
}

void bv_or(const binary_vector_t *a, const binary_vector_t *b, binary_vector_t *r) {
    if (!a || !b || !r) return;
    if (a->dimension != b->dimension || a->dimension != r->dimension) return;
    int32_t words = (a->dimension + 31) / 32;
    for (int32_t i = 0; i < words; i++)
        r->bits[i] = a->bits[i] | b->bits[i];
}

void bv_xor(const binary_vector_t *a, const binary_vector_t *b, binary_vector_t *r) {
    if (!a || !b || !r) return;
    if (a->dimension != b->dimension || a->dimension != r->dimension) return;
    int32_t words = (a->dimension + 31) / 32;
    for (int32_t i = 0; i < words; i++)
        r->bits[i] = a->bits[i] ^ b->bits[i];
}

void bv_not(const binary_vector_t *v, binary_vector_t *r) {
    if (!v || !r || v->dimension != r->dimension) return;
    int32_t words = (v->dimension + 31) / 32;
    for (int32_t i = 0; i < words; i++)
        r->bits[i] = ~v->bits[i];
    int32_t rem = v->dimension % 32;
    if (rem > 0) r->bits[words - 1] &= (1U << rem) - 1;
}

bool bv_is_zero(const binary_vector_t *v) {
    if (!v) return true;
    int32_t words = (v->dimension + 31) / 32;
    for (int32_t i = 0; i < words; i++)
        if (v->bits[i] != 0) return false;
    return true;
}

bool bv_is_ones(const binary_vector_t *v) {
    if (!v) return false;
    int32_t words = v->dimension / 32;
    for (int32_t i = 0; i < words; i++)
        if (v->bits[i] != 0xFFFFFFFFU) return false;
    int32_t rem = v->dimension % 32;
    if (rem > 0 && v->bits[words] != ((1U << rem) - 1)) return false;
    return true;
}

bool bv_equals(const binary_vector_t *a, const binary_vector_t *b) {
    if (!a || !b) return (a == b);
    if (a->dimension != b->dimension) return false;
    int32_t words = (a->dimension + 31) / 32;
    for (int32_t i = 0; i < words; i++)
        if (a->bits[i] != b->bits[i]) return false;
    return true;
}

int32_t bv_popcount32(uint32_t w) {
    return (int32_t)BV_POPCOUNT32(w);
}

bool bv_is_subset(const binary_vector_t *v, const binary_vector_t *w) {
    if (!v || !w || v->dimension != w->dimension) return false;
    int32_t words = (v->dimension + 31) / 32;
    for (int32_t i = 0; i < words; i++)
        if ((v->bits[i] & ~w->bits[i]) != 0) return false;
    return true;
}

/* ===================================================================
 * Utilities
 * =================================================================== */

double ov_get_time_seconds(void) {
    return (double)clock() / (double)CLOCKS_PER_SEC;
}

void ov_print_summary(const ov_instance_t *inst) {
    if (!inst) { printf("OV Instance: NULL\n"); return; }
    printf("=== OV Instance ===\n");
    printf("  |A| = %d, |B| = %d\n",
           inst->A ? inst->A->num_vectors : 0,
           inst->B ? inst->B->num_vectors : 0);
    printf("  dim = %d\n", inst->A ? inst->A->dimension : 0);
    printf("  Total pairs: %lld\n", (long long)ov_num_pairs(inst));
    if (inst->A) {
        density_stats_t ds = vs_density_stats(inst->A);
        printf("  A density: mean=%.3f min=%.3f max=%.3f\n",
               ds.mean_density, ds.min_density, ds.max_density);
    }
}

bool ov_validate(const ov_instance_t *inst) {
    if (!inst) return false;
    if (!inst->A || !inst->B) return false;
    if (inst->A->dimension != inst->B->dimension) return false;
    if (inst->A->num_vectors <= 0 || inst->B->num_vectors <= 0) return false;
    for (int32_t i = 0; i < inst->A->num_vectors; i++)
        if (inst->A->vectors[i].dimension != inst->A->dimension) return false;
    for (int32_t i = 0; i < inst->B->num_vectors; i++)
        if (inst->B->vectors[i].dimension != inst->B->dimension) return false;
    return true;
}

void ov_report_result(const ov_result_t *r) {
    if (!r) return;
    printf("OV Result: found=%s", r->found ? "YES" : "NO");
    if (r->found) printf(" a[%d] * b[%d] = 0", r->a_index, r->b_index);
    printf(" ops=%lld time=%.6fs\n",
           (long long)r->ops, r->time_s);
}
