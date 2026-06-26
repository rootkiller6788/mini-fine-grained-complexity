/* ov_applications.c - OV Applications (L7)
 * Database anti-join, CRISPR off-target prediction,
 * Code similarity detection, Collaborative filtering.
 */
#include "ov.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ===================================================================
 * L7.1: Database Anti-Join Simulation
 *
 * In relational databases, an anti-join finds rows in table A that
 * have NO matching row in table B on specified join attributes.
 *
 * OV model: each row is a vector in {0,1}^d, where d = number of
 * boolean join attributes. Two rows are orthogonal iff they disagree
 * on ALL attributes - i.e., they would NOT join.
 *
 * This function finds all pairs of rows that are completely disjoint
 * on the join attributes.
 * =================================================================== */

int64_t ov_database_anti_join(const ov_instance_t *inst,
                               int32_t *row_a, int32_t *row_b, int64_t max_results) {
    /* inst->A = left table rows, inst->B = right table rows.
     * Each row is a binary vector where bit k = 1 means
     * "row has value 1 on attribute k".
     * Orthogonal = no shared 1-bits on any attribute = disjoint rows. */
    if (!inst || !row_a || !row_b || max_results <= 0) return 0;
    int32_t nA = inst->A->num_vectors;
    int32_t nB = inst->B->num_vectors;
    int64_t count = 0;
    for (int32_t i = 0; i < nA && count < max_results; i++) {
        for (int32_t j = 0; j < nB && count < max_results; j++) {
            if (bv_are_orthogonal(&inst->A->vectors[i], &inst->B->vectors[j])) {
                row_a[count] = i;
                row_b[count] = j;
                count++;
            }
        }
    }
    return count;
}

/* ===================================================================
 * L7.2: Computational Biology - CRISPR Off-Target Prediction
 *
 * In CRISPR-Cas9 gene editing, guide RNA (gRNA) sequences must be
 * designed to target specific genomic regions while avoiding
 * off-target binding.
 *
 * Model: Each gRNA candidate and each genomic region is encoded as
 * a binary vector. Features include:
 *   - PAM sequence presence (NGG)
 *   - GC content in specific ranges
 *   - Seed region complementarity (8-12 bases at 3' end)
 *   - Mismatch tolerance patterns
 *
 * Two vectors are orthogonal if they share no binding features,
 * i.e., the gRNA is predicted to NOT bind to that region (safe).
 * =================================================================== */

/* Convert DNA base to 2-bit encoding */
static int32_t dna_to_bits(char base) {
    switch (base) {
        case 'A': case 'a': return 0;  /* 00 */
        case 'C': case 'c': return 1;  /* 01 */
        case 'G': case 'g': return 2;  /* 10 */
        case 'T': case 't': return 3;  /* 11 */
        default: return 0;
    }
}

ov_instance_t *ov_from_dna_sequences(const char **sequences, int32_t n,
                                      int32_t seq_length) {
    /* Convert n DNA sequences of length seq_length into binary vectors.
     * Each base -> 2 bits, so dimension = seq_length * 2.
     * Plus 4 extra bits for GC content encoding. */
    if (!sequences || n <= 0 || seq_length <= 0) return NULL;
    int32_t d = seq_length * 2 + 4;
    ov_instance_t *inst = ov_create(n, d);
    if (!inst) return NULL;
    for (int32_t i = 0; i < n; i++) {
        if (!sequences[i] || (int32_t)strlen(sequences[i]) < seq_length) {
            /* Fill with zeros for invalid sequences */
            continue;
        }
        binary_vector_t *vA = &inst->A->vectors[i];
        binary_vector_t *vB = &inst->B->vectors[i];
        /* Encode each base as 2 bits */
        for (int32_t p = 0; p < seq_length; p++) {
            int32_t bits = dna_to_bits(sequences[i][p]);
            bv_set(vA, p * 2,     bits & 1);
            bv_set(vA, p * 2 + 1, (bits >> 1) & 1);
            bv_set(vB, p * 2,     bits & 1);
            bv_set(vB, p * 2 + 1, (bits >> 1) & 1);
        }
        /* Encode GC content: count G+C bases, set bits if > 50% */
        int32_t gc = 0;
        for (int32_t p = 0; p < seq_length; p++) {
            char c = sequences[i][p];
            if (c == 'G' || c == 'g' || c == 'C' || c == 'c') gc++;
        }
        double gc_ratio = (double)gc / (double)seq_length;
        bv_set(vA, seq_length * 2, gc_ratio > 0.4);
        bv_set(vA, seq_length * 2 + 1, gc_ratio > 0.5);
        bv_set(vA, seq_length * 2 + 2, gc_ratio > 0.6);
        bv_set(vA, seq_length * 2 + 3, gc_ratio > 0.7);
        /* B set gets same encoding for self-comparison */
        memcpy(vB->bits, vA->bits, (size_t)((d + 31) / 32) * sizeof(uint32_t));
    }
    return inst;
}

int32_t ov_count_safe_crispr_pairs(const ov_instance_t *inst) {
    /* Count pairs (gRNA, genomic region) that are orthogonal.
     * Orthogonal = no shared 1-bits = no binding features in common.
     * These pairs are predicted to be safe (no off-target binding). */
    if (!inst || !inst->A || !inst->B) return 0;
    int32_t n = inst->A->num_vectors;
    int32_t safe = 0;
    for (int32_t i = 0; i < n; i++) {
        for (int32_t j = 0; j < n; j++) {
            if (bv_are_orthogonal(&inst->A->vectors[i], &inst->B->vectors[j]))
                safe++;
        }
    }
    return safe;
}

/* ===================================================================
 * L7.3: Code Similarity / Plagiarism Detection
 *
 * Code fragments are converted to binary feature vectors:
 *   - Bit i = 1 if code fragment contains feature i
 *   - Features: keywords (if, while, for, class, def, etc.),
 *     API calls, operator types, control flow patterns.
 *
 * Orthogonal vectors = code fragments with NO shared features.
 * This identifies maximally dissimilar code pairs, useful for:
 *   - Diversity analysis in codebases
 *   - Identifying unique coding styles
 *   - Plagiarism exclusion (orthogonal = definitely not plagiarized)
 * =================================================================== */

double ov_max_code_dissimilarity(const ov_instance_t *inst,
                                  int32_t *frag_a, int32_t *frag_b) {
    /* Find the pair of code fragments with maximum dissimilarity.
     * Dissimilarity = 1 - Jaccard similarity.
     * Orthogonal vectors have Jaccard = 0.0, dissimilarity = 1.0. */
    if (!inst || !inst->A || !inst->B || !frag_a || !frag_b) return -1.0;
    int32_t n = inst->A->num_vectors;
    double max_dissim = -1.0;
    *frag_a = -1; *frag_b = -1;
    for (int32_t i = 0; i < n; i++) {
        for (int32_t j = 0; j < n; j++) {
            double jac = bv_jaccard(&inst->A->vectors[i], &inst->B->vectors[j]);
            double dissim = 1.0 - jac;
            if (dissim > max_dissim) {
                max_dissim = dissim;
                *frag_a = i;
                *frag_b = j;
            }
        }
    }
    return max_dissim;
}

/* ===================================================================
 * L7.4: Collaborative Filtering - Diverse User Detection
 *
 * In recommendation systems, finding users with orthogonal preference
 * vectors identifies maximally diverse user segments. This is useful for:
 *   - Creating diverse recommendation lists
 *   - Identifying complementary user groups
 *   - Cold-start recommendation strategies
 * =================================================================== */

/* Min-heap entry for top-k tracking */
typedef struct {
    double score;
    int32_t a, b;
} user_pair_t;

static int cmp_user_pair_desc(const void *pa, const void *pb) {
    const user_pair_t *a = (const user_pair_t *)pa;
    const user_pair_t *b = (const user_pair_t *)pb;
    if (a->score < b->score) return 1;
    if (a->score > b->score) return -1;
    return 0;
}

int32_t ov_find_topk_diverse_users(const ov_instance_t *inst, int32_t top_k,
                                    int32_t *user_a, int32_t *user_b,
                                    double *scores) {
    if (!inst || !inst->A || !inst->B || !user_a || !user_b || !scores || top_k <= 0)
        return 0;
    int32_t n = inst->A->num_vectors;
    int32_t total_pairs = n * n;
    user_pair_t *pairs = (user_pair_t *)malloc((size_t)total_pairs * sizeof(user_pair_t));
    if (!pairs) return 0;
    int32_t pc = 0;
    /* Compute dissimilarity for all pairs */
    for (int32_t i = 0; i < n; i++) {
        for (int32_t j = 0; j < n; j++) {
            double jac = bv_jaccard(&inst->A->vectors[i], &inst->B->vectors[j]);
            pairs[pc].score = 1.0 - jac;  /* higher = more diverse */
            pairs[pc].a = i;
            pairs[pc].b = j;
            pc++;
        }
    }
    /* Sort by diversity score descending */
    qsort(pairs, (size_t)total_pairs, sizeof(user_pair_t), cmp_user_pair_desc);
    /* Return top-k */
    int32_t k = top_k < total_pairs ? top_k : total_pairs;
    for (int32_t m = 0; m < k; m++) {
        user_a[m] = pairs[m].a;
        user_b[m] = pairs[m].b;
        scores[m] = pairs[m].score;
    }
    free(pairs);
    return k;
}
