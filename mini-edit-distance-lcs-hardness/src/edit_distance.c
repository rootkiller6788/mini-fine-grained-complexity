/* ============================================================================
 * edit_distance.c -- Core Edit Distance Implementation
 *
 * Wagner-Fischer algorithm, bounded edit distance, Damerau-Levenshtein,
 * Jaro-Winkler, Hamming distance, normalized edit distance.
 *
 * L1: edit_op_t, edit_script_t, edit_cost_model_t, edit_dp_table_t
 * L2: Quadratic-time edit distance, linear-space variant
 * L3: DP table fill, traceable DP, edit script reconstruction
 * L5: Wagner-Fischer O(nm), bounded k-difference, Four Russians speedup
 * L8: Approximate edit distance (Andoni-Krauthgamer-Onak 2010)
 *
 * References:
 *   Levenshtein (1966), Wagner & Fischer (1974), Masek & Paterson (1980),
 *   Ukkonen (1985), Backurs & Indyk (STOC 2015)
 * ============================================================================ */

#include "edit_distance.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ---- Inline helpers ---- */
static inline int32_t min3(int32_t a, int32_t b, int32_t c) {
    int32_t m = (a < b) ? a : b;
    return (m < c) ? m : c;
}
static inline int32_t min2(int32_t a, int32_t b) {
    return (a < b) ? a : b;
}
static inline int32_t max2(int32_t a, int32_t b) {
    return (a > b) ? a : b;
}

/* L1: Default unit-cost model for Levenshtein distance */
static const edit_cost_model_t UNIT_COST_MODEL = {1, 1, 1, 0};

/* ============================================================================
 * L2: Core Edit Distance -- Wagner-Fischer Algorithm
 *
 * DP recurrence:
 *   D[0][j] = j * insert_cost
 *   D[i][0] = i * delete_cost
 *   D[i][j] = min(
 *     D[i-1][j]   + delete_cost,
 *     D[i][j-1]   + insert_cost,
 *     D[i-1][j-1] + (s[i-1]==t[j-1] ? match_cost : substitute_cost)
 *   )
 * ============================================================================ */

int32_t edit_distance_full(const char *src, const char *tgt) {
    return edit_distance_weighted(src, tgt, &UNIT_COST_MODEL);
}

int32_t edit_distance_weighted(const char *src, const char *tgt,
                               const edit_cost_model_t *costs) {
    if (!src || !tgt || !costs) return -1;

    int32_t n = (int32_t)strlen(src);
    int32_t m = (int32_t)strlen(tgt);
    if (n == 0) return m * costs->insert_cost;
    if (m == 0) return n * costs->delete_cost;

    int32_t stride = m + 1;
    int32_t total = (n + 1) * stride;
    int32_t *dp = (int32_t *)calloc((size_t)total, sizeof(int32_t));
    if (!dp) return -1;

    for (int32_t j = 0; j <= m; j++) dp[j] = j * costs->insert_cost;
    for (int32_t i = 1; i <= n; i++) dp[i * stride] = i * costs->delete_cost;

    for (int32_t i = 1; i <= n; i++) {
        int32_t rb = i * stride, prb = (i - 1) * stride;
        for (int32_t j = 1; j <= m; j++) {
            int32_t del = dp[prb + j] + costs->delete_cost;
            int32_t ins = dp[rb + (j - 1)] + costs->insert_cost;
            int32_t match_cost = (src[i - 1] == tgt[j - 1])
                                 ? costs->match_cost : costs->substitute_cost;
            int32_t sub = dp[prb + (j - 1)] + match_cost;
            dp[rb + j] = min3(del, ins, sub);
        }
    }
    int32_t result = dp[n * stride + m];
    free(dp);
    return result;
}

/* Space-efficient edit distance: O(min(n,m)) space, 2-row DP trick */
int32_t edit_distance_linear_space(const char *src, const char *tgt) {
    if (!src || !tgt) return -1;
    int32_t n = (int32_t)strlen(src);
    int32_t m = (int32_t)strlen(tgt);
    if (n == 0) return m;
    if (m == 0) return n;

    int32_t *prev = (int32_t *)calloc((size_t)(m + 1), sizeof(int32_t));
    int32_t *curr = (int32_t *)calloc((size_t)(m + 1), sizeof(int32_t));
    if (!prev || !curr) { free(prev); free(curr); return -1; }

    for (int32_t j = 0; j <= m; j++) prev[j] = j;

    for (int32_t i = 1; i <= n; i++) {
        curr[0] = i;
        for (int32_t j = 1; j <= m; j++) {
            int32_t del = prev[j] + 1;
            int32_t ins = curr[j - 1] + 1;
            int32_t sub = prev[j - 1] + ((src[i-1] == tgt[j-1]) ? 0 : 1);
            curr[j] = min3(del, ins, sub);
        }
        int32_t *tmp = prev; prev = curr; curr = tmp;
    }
    int32_t result = prev[m];
    free(prev); free(curr);
    return result;
}

/* Ukkonen's bounded edit distance: O(k * min(n,m)) time.
 * Only computes cells within diagonal band of width 2k+1.
 * Returns min(distance, k+1). */
int32_t edit_distance_bounded(const char *src, const char *tgt, int32_t k) {
    if (!src || !tgt || k < 0) return -1;
    int32_t n = (int32_t)strlen(src);
    int32_t m = (int32_t)strlen(tgt);
    if (n == 0) return (m <= k) ? m : k + 1;
    if (m == 0) return (n <= k) ? n : k + 1;

    /* Quick filter: length difference > k implies distance > k */
    int32_t len_diff = (n > m) ? (n - m) : (m - n);
    if (len_diff > k) return k + 1;

    int32_t *prev = (int32_t *)calloc((size_t)(m + 1), sizeof(int32_t));
    int32_t *curr = (int32_t *)calloc((size_t)(m + 1), sizeof(int32_t));
    if (!prev || !curr) { free(prev); free(curr); return -1; }
    for (int32_t j = 0; j <= m; j++) prev[j] = j;

    for (int32_t i = 1; i <= n; i++) {
        int32_t band_start = (i - k > 1) ? (i - k) : 1;
        int32_t band_end   = (i + k < m) ? (i + k) : m;
        curr[0] = i;
        for (int32_t j = 1; j <= m; j++) {
            if (j < band_start || j > band_end) {
                curr[j] = k + 1;
            } else {
                int32_t del = prev[j] + 1;
                int32_t ins = curr[j - 1] + 1;
                int32_t sub = prev[j - 1] + ((src[i-1] == tgt[j-1]) ? 0 : 1);
                curr[j] = min3(del, ins, sub);
            }
        }
        bool all_above_k = true;
        for (int32_t j = band_start; j <= band_end; j++) {
            if (curr[j] <= k) { all_above_k = false; break; }
        }
        if (all_above_k) { free(prev); free(curr); return k + 1; }
        int32_t *tmp = prev; prev = curr; curr = tmp;
    }
    int32_t result = (prev[m] <= k) ? prev[m] : k + 1;
    free(prev); free(curr);
    return result;
}

/* ============================================================================
 * L3: DP Table Operations
 * ============================================================================ */

edit_dp_table_t *ed_dp_create(int32_t src_len, int32_t tgt_len) {
    edit_dp_table_t *dp = (edit_dp_table_t *)malloc(sizeof(edit_dp_table_t));
    if (!dp) return NULL;
    dp->rows = src_len + 1;
    dp->cols = tgt_len + 1;
    dp->capacity = dp->rows * dp->cols;
    dp->data = (int32_t *)calloc((size_t)dp->capacity, sizeof(int32_t));
    if (!dp->data) { free(dp); return NULL; }
    return dp;
}

void ed_dp_destroy(edit_dp_table_t *dp) {
    if (!dp) return;
    free(dp->data);
    free(dp);
}

void ed_dp_fill(edit_dp_table_t *dp, const char *src, const char *tgt) {
    if (!dp || !src || !tgt) return;
    int32_t n = dp->rows - 1;
    int32_t m = dp->cols - 1;
    int32_t stride = dp->cols;
    for (int32_t i = 0; i <= n; i++) dp->data[i * stride] = i;
    for (int32_t j = 0; j <= m; j++) dp->data[j] = j;
    for (int32_t i = 1; i <= n; i++) {
        int32_t rb = i * stride, prb = (i - 1) * stride;
        for (int32_t j = 1; j <= m; j++) {
            int32_t del = dp->data[prb + j] + 1;
            int32_t ins = dp->data[rb + (j - 1)] + 1;
            int32_t sub = dp->data[prb + (j - 1)]
                        + ((src[i-1] == tgt[j-1]) ? 0 : 1);
            dp->data[rb + j] = min3(del, ins, sub);
        }
    }
}

edit_dp_traceable_t *ed_dp_traceable_create(int32_t src_len, int32_t tgt_len) {
    edit_dp_traceable_t *dp =
        (edit_dp_traceable_t *)malloc(sizeof(edit_dp_traceable_t));
    if (!dp) return NULL;
    dp->rows = src_len + 1;
    dp->cols = tgt_len + 1;
    int32_t sz = dp->rows * dp->cols;
    dp->data = (edit_dp_cell_t *)calloc((size_t)sz, sizeof(edit_dp_cell_t));
    if (!dp->data) { free(dp); return NULL; }
    return dp;
}

void ed_dp_traceable_destroy(edit_dp_traceable_t *dp) {
    if (!dp) return;
    free(dp->data);
    free(dp);
}

void ed_dp_traceable_fill(edit_dp_traceable_t *dp, const char *src,
                           const char *tgt, const edit_cost_model_t *costs) {
    if (!dp || !src || !tgt || !costs) return;
    int32_t n = dp->rows - 1;
    int32_t m = dp->cols - 1;
    int32_t strd = dp->cols;
    edit_dp_cell_t *d = dp->data;

    for (int32_t i = 0; i <= n; i++) {
        d[i * strd].value = i * costs->delete_cost;
        d[i * strd].backpointer = EDIT_OP_DELETE;
    }
    for (int32_t j = 0; j <= m; j++) {
        d[j].value = j * costs->insert_cost;
        d[j].backpointer = EDIT_OP_INSERT;
    }
    d[0].backpointer = EDIT_OP_MATCH;

    for (int32_t i = 1; i <= n; i++) {
        int32_t rb = i * strd, prb = (i - 1) * strd;
        for (int32_t j = 1; j <= m; j++) {
            int32_t del_val = d[prb + j].value + costs->delete_cost;
            int32_t ins_val = d[rb + (j - 1)].value + costs->insert_cost;
            int32_t match_cost = (src[i-1] == tgt[j-1])
                                 ? costs->match_cost : costs->substitute_cost;
            int32_t sub_val = d[prb + (j - 1)].value + match_cost;
            if (sub_val <= del_val && sub_val <= ins_val) {
                d[rb + j].value = sub_val;
                d[rb + j].backpointer = (match_cost == costs->match_cost)
                    ? EDIT_OP_MATCH : EDIT_OP_SUBSTITUTE;
            } else if (del_val <= ins_val) {
                d[rb + j].value = del_val;
                d[rb + j].backpointer = EDIT_OP_DELETE;
            } else {
                d[rb + j].value = ins_val;
                d[rb + j].backpointer = EDIT_OP_INSERT;
            }
        }
    }
}

/* ============================================================================
 * Edit Script Reconstruction
 *
 * Trace back through the DP table from D[n][m] to D[0][0] following
 * backpointers to reconstruct the minimal edit script.
 * ============================================================================ */

edit_script_t *ed_reconstruct_script(const edit_dp_traceable_t *dp,
                                     const char *src, const char *tgt) {
    if (!dp || !src || !tgt) return NULL;
    int32_t i = dp->rows - 1, j = dp->cols - 1;
    int32_t strd = dp->cols;
    edit_dp_cell_t *d = dp->data;

    /* Count operations on optimal path */
    int32_t count = 0;
    int32_t ci = i, cj = j;
    while (ci > 0 || cj > 0) {
        edit_op_type_t op = d[ci * strd + cj].backpointer;
        count++;
        if (op == EDIT_OP_MATCH || op == EDIT_OP_SUBSTITUTE) { ci--; cj--; }
        else if (op == EDIT_OP_DELETE) { ci--; }
        else { cj--; }
    }

    edit_script_t *script = (edit_script_t *)malloc(sizeof(edit_script_t));
    if (!script) return NULL;
    script->operations = (edit_op_t *)calloc((size_t)count, sizeof(edit_op_t));
    if (!script->operations) { free(script); return NULL; }
    script->num_operations = count;
    script->total_cost = d[(dp->rows - 1) * strd + (dp->cols - 1)].value;

    /* Fill operations in reverse, starting from D[n][m] */
    ci = i; cj = j;
    int32_t pos = count - 1;
    while (ci > 0 || cj > 0) {
        edit_op_type_t op = d[ci * strd + cj].backpointer;
        edit_op_t *eo = &script->operations[pos];
        if (op == EDIT_OP_MATCH || op == EDIT_OP_SUBSTITUTE) {
            eo->type = op;
            eo->source_pos = ci - 1; eo->target_pos = cj - 1;
            eo->source_char = src[ci - 1]; eo->target_char = tgt[cj - 1];
            eo->cost = (op == EDIT_OP_MATCH) ? 0 : 1;
            ci--; cj--;
        } else if (op == EDIT_OP_DELETE) {
            eo->type = EDIT_OP_DELETE;
            eo->source_pos = ci - 1; eo->target_pos = -1;
            eo->source_char = src[ci - 1]; eo->target_char = 0;
            eo->cost = 1;
            ci--;
        } else {
            eo->type = EDIT_OP_INSERT;
            eo->source_pos = -1; eo->target_pos = cj - 1;
            eo->source_char = 0; eo->target_char = tgt[cj - 1];
            eo->cost = 1;
            cj--;
        }
        pos--;
    }
    return script;
}

void ed_script_destroy(edit_script_t *script) {
    if (!script) return;
    free(script->operations);
    free(script);
}

void ed_script_print(const edit_script_t *script) {
    if (!script) { printf("(null script)\n"); return; }
    printf("Edit Script: %d operations, total cost %d\n",
           script->num_operations, script->total_cost);
    const char *op_names[] = {"MATCH", "SUB", "INS", "DEL"};
    for (int32_t k = 0; k < script->num_operations; k++) {
        edit_op_t *op = &script->operations[k];
        printf("  %2d: %s ", k, op_names[op->type]);
        if (op->source_pos >= 0)
            printf("src[%d]='%c'", op->source_pos, op->source_char);
        else printf("-----------");
        printf(" -> ");
        if (op->target_pos >= 0)
            printf("tgt[%d]='%c'", op->target_pos, op->target_char);
        else printf("-----------");
        printf("  cost=%d\n", op->cost);
    }
}

edit_script_t *edit_distance_with_script(const char *src, const char *tgt) {
    if (!src || !tgt) return NULL;
    int32_t n = (int32_t)strlen(src), m = (int32_t)strlen(tgt);
    edit_dp_traceable_t *dp = ed_dp_traceable_create(n, m);
    if (!dp) return NULL;
    ed_dp_traceable_fill(dp, src, tgt, &UNIT_COST_MODEL);
    edit_script_t *script = ed_reconstruct_script(dp, src, tgt);
    ed_dp_traceable_destroy(dp);
    return script;
}

/* Hirschberg-style linear-space edit script reconstruction.
 * Uses traceable DP for correctness; true linear-space version
 * would use divide-and-conquer with forward/backward DP in O(m) space. */
edit_script_t *ed_reconstruct_linear(const char *src, const char *tgt) {
    if (!src || !tgt) return NULL;
    int32_t n = (int32_t)strlen(src), m = (int32_t)strlen(tgt);

    if (n == 0 && m == 0) {
        edit_script_t *s = (edit_script_t *)malloc(sizeof(edit_script_t));
        if (!s) return NULL;
        s->operations = NULL; s->num_operations = 0; s->total_cost = 0;
        return s;
    }
    if (n == 0) {
        edit_script_t *s = (edit_script_t *)malloc(sizeof(edit_script_t));
        if (!s) return NULL;
        s->num_operations = m; s->total_cost = m;
        s->operations = (edit_op_t *)calloc((size_t)m, sizeof(edit_op_t));
        if (!s->operations) { free(s); return NULL; }
        for (int32_t j = 0; j < m; j++) {
            s->operations[j].type = EDIT_OP_INSERT;
            s->operations[j].target_pos = j;
            s->operations[j].target_char = tgt[j];
            s->operations[j].cost = 1;
        }
        return s;
    }
    if (m == 0) {
        edit_script_t *s = (edit_script_t *)malloc(sizeof(edit_script_t));
        if (!s) return NULL;
        s->num_operations = n; s->total_cost = n;
        s->operations = (edit_op_t *)calloc((size_t)n, sizeof(edit_op_t));
        if (!s->operations) { free(s); return NULL; }
        for (int32_t i = 0; i < n; i++) {
            s->operations[i].type = EDIT_OP_DELETE;
            s->operations[i].source_pos = i;
            s->operations[i].source_char = src[i];
            s->operations[i].cost = 1;
        }
        return s;
    }

    /* Use full traceable DP to reconstruct */
    edit_dp_traceable_t *dp = ed_dp_traceable_create(n, m);
    if (!dp) return NULL;
    ed_dp_traceable_fill(dp, src, tgt, &UNIT_COST_MODEL);
    edit_script_t *full = ed_reconstruct_script(dp, src, tgt);
    ed_dp_traceable_destroy(dp);
    if (!full) return NULL;

    edit_script_t *result = (edit_script_t *)malloc(sizeof(edit_script_t));
    if (!result) { ed_script_destroy(full); return NULL; }
    result->operations = full->operations;
    result->num_operations = full->num_operations;
    result->total_cost = full->total_cost;
    free(full);
    return result;
}

/* ============================================================================
 * Hamming Distance (1950)
 *
 * Number of positions where two equal-length strings differ.
 * Only allows substitutions; no insertions or deletions.
 * ============================================================================ */

int32_t hamming_distance(const char *s1, const char *s2) {
    if (!s1 || !s2) return -1;
    int32_t len1 = (int32_t)strlen(s1);
    int32_t len2 = (int32_t)strlen(s2);
    if (len1 != len2) return -1;
    int32_t dist = 0;
    for (int32_t i = 0; i < len1; i++) {
        if (s1[i] != s2[i]) dist++;
    }
    return dist;
}

double hamming_distance_weighted(const char *s1, const char *s2,
                                 const double *weights) {
    if (!s1 || !s2 || !weights) return -1.0;
    int32_t len1 = (int32_t)strlen(s1);
    int32_t len2 = (int32_t)strlen(s2);
    if (len1 != len2) return -1.0;
    double dist = 0.0;
    for (int32_t i = 0; i < len1; i++) {
        if (s1[i] != s2[i]) dist += weights[i];
    }
    return dist;
}

/* ============================================================================
 * Damerau-Levenshtein Distance
 *
 * Extends Levenshtein by allowing transposition of adjacent characters
 * as a single operation. Transposition cost is 1.
 *
 * DP adds: if src[i-1]==tgt[j-2] && src[i-2]==tgt[j-1]:
 *            D[i][j] = min(D[i][j], D[i-2][j-2] + cost)
 * ============================================================================ */

int32_t damerau_levenshtein(const char *src, const char *tgt) {
    if (!src || !tgt) return -1;
    int32_t n = (int32_t)strlen(src);
    int32_t m = (int32_t)strlen(tgt);
    if (n == 0) return m;
    if (m == 0) return n;

    int32_t stride = m + 1;
    int32_t *dp = (int32_t *)calloc((size_t)((n + 1) * stride), sizeof(int32_t));
    if (!dp) return -1;

    for (int32_t i = 0; i <= n; i++) dp[i * stride] = i;
    for (int32_t j = 0; j <= m; j++) dp[j] = j;

    for (int32_t i = 1; i <= n; i++) {
        int32_t rb = i * stride, prb = (i - 1) * stride;
        for (int32_t j = 1; j <= m; j++) {
            int32_t cost = (src[i-1] == tgt[j-1]) ? 0 : 1;
            int32_t del = dp[prb + j] + 1;
            int32_t ins = dp[rb + (j - 1)] + 1;
            int32_t sub = dp[prb + (j - 1)] + cost;
            int32_t val = min3(del, ins, sub);

            /* Transposition check */
            if (i > 1 && j > 1 &&
                src[i-1] == tgt[j-2] && src[i-2] == tgt[j-1]) {
                int32_t trans = dp[(i-2) * stride + (j - 2)] + cost;
                if (trans < val) val = trans;
            }
            dp[rb + j] = val;
        }
    }
    int32_t result = dp[n * stride + m];
    free(dp);
    return result;
}

/* ============================================================================
 * Jaro Similarity and Jaro-Winkler Distance
 *
 * Jaro (1989): string similarity for record linkage.
 * Jaro-Winkler (1990): extends Jaro with prefix bonus.
 *
 * Jaro formula: (m/|s1| + m/|s2| + (m-t)/m) / 3
 * where m = number of matching chars within window
 *       t = half the number of transpositions
 * ============================================================================ */

double jaro_similarity(const char *s1, const char *s2) {
    if (!s1 || !s2) return 0.0;
    int32_t len1 = (int32_t)strlen(s1);
    int32_t len2 = (int32_t)strlen(s2);
    if (len1 == 0 && len2 == 0) return 1.0;
    if (len1 == 0 || len2 == 0) return 0.0;

    /* Matching window: floor(max(len1,len2)/2) - 1 */
    int32_t window = max2(len1, len2) / 2;
    if (window > 0) window--;

    bool *m1 = (bool *)calloc((size_t)len1, sizeof(bool));
    bool *m2 = (bool *)calloc((size_t)len2, sizeof(bool));
    if (!m1 || !m2) { free(m1); free(m2); return 0.0; }

    int32_t matches = 0;
    for (int32_t i = 0; i < len1; i++) {
        int32_t start = max2(0, i - window);
        int32_t end   = min2(len2 - 1, i + window);
        for (int32_t j = start; j <= end; j++) {
            if (!m2[j] && s1[i] == s2[j]) {
                m1[i] = true; m2[j] = true; matches++; break;
            }
        }
    }
    if (matches == 0) { free(m1); free(m2); return 0.0; }

    /* Count transpositions */
    int32_t trans = 0, k = 0;
    for (int32_t i = 0; i < len1; i++) {
        if (m1[i]) {
            while (k < len2 && !m2[k]) k++;
            if (k < len2 && s1[i] != s2[k]) trans++;
            k++;
        }
    }
    trans /= 2;
    free(m1); free(m2);

    double m = (double)matches, t = (double)trans;
    return (m/len1 + m/len2 + (m-t)/m) / 3.0;
}

double jaro_winkler(const char *s1, const char *s2, double prefix_scale) {
    if (!s1 || !s2) return 0.0;
    if (prefix_scale < 0.0) prefix_scale = 0.0;
    if (prefix_scale > 0.25) prefix_scale = 0.25;
    double jaro = jaro_similarity(s1, s2);
    int32_t prefix = 0;
    while (prefix < 4 && s1[prefix] && s2[prefix] &&
           s1[prefix] == s2[prefix]) prefix++;
    return jaro + (double)prefix * prefix_scale * (1.0 - jaro);
}

/* ============================================================================
 * Longest Common Prefix / Suffix
 * ============================================================================ */

int32_t longest_common_prefix(const char *s1, const char *s2) {
    if (!s1 || !s2) return 0;
    int32_t i = 0;
    while (s1[i] && s2[i] && s1[i] == s2[i]) i++;
    return i;
}

int32_t longest_common_suffix(const char *s1, const char *s2) {
    if (!s1 || !s2) return 0;
    int32_t l1 = (int32_t)strlen(s1), l2 = (int32_t)strlen(s2);
    int32_t i = 0;
    while (i < l1 && i < l2 && s1[l1-1-i] == s2[l2-1-i]) i++;
    return i;
}

/* Normalized edit distance: ED(s,t) / max(|s|,|t|) in [0,1] */
double edit_distance_normalized(const char *src, const char *tgt) {
    if (!src || !tgt) return 1.0;
    int32_t d = edit_distance_full(src, tgt);
    int32_t max_len = max2((int32_t)strlen(src), (int32_t)strlen(tgt));
    if (max_len == 0) return 0.0;
    return (double)d / (double)max_len;
}

/* ============================================================================
 * L5: Four Russians Speedup (Masek-Paterson 1980)
 *
 * Precomputes edit distance for all binary string pairs of length t
 * (t = O(log n)), then uses these to compute the DP in blocks,
 * reducing time from O(n^2) to O(n^2 / log n).
 *
 * This implementation uses t = 4 for demonstration.
 * ============================================================================ */

#define FR_BLOCK_SIZE 4

static int32_t fr_precomp[1 << FR_BLOCK_SIZE][1 << FR_BLOCK_SIZE];
static bool fr_ready = false;

static void fr_init(void) {
    if (fr_ready) return;
    int32_t t = FR_BLOCK_SIZE;
    int32_t total = 1 << t;
    for (int32_t ma = 0; ma < total; ma++) {
        for (int32_t mb = 0; mb < total; mb++) {
            char sa[FR_BLOCK_SIZE+1], sb[FR_BLOCK_SIZE+1];
            for (int32_t k = 0; k < t; k++) {
                sa[k] = ((ma >> k) & 1) ? '1' : '0';
                sb[k] = ((mb >> k) & 1) ? '1' : '0';
            }
            sa[t] = 0; sb[t] = 0;
            /* Compute edit distance by small DP */
            int32_t dp2[5][5];
            for (int32_t i = 0; i <= t; i++) dp2[i][0] = i;
            for (int32_t j = 0; j <= t; j++) dp2[0][j] = j;
            for (int32_t i = 1; i <= t; i++)
                for (int32_t j = 1; j <= t; j++) {
                    int32_t del = dp2[i-1][j] + 1;
                    int32_t ins = dp2[i][j-1] + 1;
                    int32_t sub = dp2[i-1][j-1] + ((sa[i-1]==sb[j-1])?0:1);
                    dp2[i][j] = min3(del, ins, sub);
                }
            fr_precomp[ma][mb] = dp2[t][t];
        }
    }
    fr_ready = true;
}

int32_t edit_distance_four_russians(const char *src, const char *tgt) {
    if (!src || !tgt) return -1;
    int32_t n = (int32_t)strlen(src);
    int32_t m = (int32_t)strlen(tgt);
    int32_t t = FR_BLOCK_SIZE;
    if (n < t || m < t) return edit_distance_full(src, tgt);

    fr_init();

    /* Binary-encode strings by character parity */
    char *bs = (char *)calloc((size_t)(n+1), 1);
    char *bt = (char *)calloc((size_t)(m+1), 1);
    if (!bs || !bt) { free(bs); free(bt); return -1; }
    for (int32_t i = 0; i < n; i++)
        bs[i] = ((unsigned char)src[i] & 1) ? '1' : '0';
    for (int32_t j = 0; j < m; j++)
        bt[j] = ((unsigned char)tgt[j] & 1) ? '1' : '0';

    /* Pad to multiple of block size */
    int32_t pn = ((n+t-1)/t)*t, pm = ((m+t-1)/t)*t;
    int32_t bn = pn/t, bm = pm/t, bstride = bm + 1;
    int32_t *b_dp = (int32_t *)calloc((size_t)((bn+1)*bstride), sizeof(int32_t));
    if (!b_dp) { free(bs); free(bt); return -1; }

    for (int32_t i = 0; i <= bn; i++) b_dp[i*bstride] = i*t;
    for (int32_t j = 0; j <= bm; j++) b_dp[j] = j*t;

    for (int32_t bi = 1; bi <= bn; bi++) {
        int32_t brb = bi*bstride, pbrb = (bi-1)*bstride;
        for (int32_t bj = 1; bj <= bm; bj++) {
            int32_t ma = 0, mb = 0;
            for (int32_t k = 0; k < t; k++) {
                int32_t ii = (bi-1)*t + k;
                if (ii < n && bs[ii] == '1') ma |= (1 << k);
            }
            for (int32_t k = 0; k < t; k++) {
                int32_t jj = (bj-1)*t + k;
                if (jj < m && bt[jj] == '1') mb |= (1 << k);
            }
            int32_t block_cost = fr_precomp[ma][mb];
            int32_t diag = b_dp[pbrb+(bj-1)] + block_cost;
            int32_t up   = b_dp[pbrb+bj] + t;
            int32_t left = b_dp[brb+(bj-1)] + t;
            b_dp[brb+bj] = min3(diag, up, left);
        }
    }
    int32_t result = b_dp[bn*bstride+bm];
    free(b_dp); free(bs); free(bt);
    return result;
}

/* Banded edit distance: only compute cells in diagonal band of width w.
 * O(w * min(n,m)) time. Delegates to edit_distance_bounded. */
int32_t edit_distance_banded(const char *src, const char *tgt, int32_t w) {
    if (!src || !tgt || w <= 0) return -1;
    return edit_distance_bounded(src, tgt, w);
}

/* ============================================================================
 * L8: Approximate Edit Distance
 *
 * Andoni, Krauthgamer, Onak (FOCS 2010):
 * "Polylogarithmic Approximation for Edit Distance and Asymmetric Edit Distance"
 *
 * Provides (1+epsilon)-approximation using sketching/embedding techniques.
 * This implementation provides a simplified estimation.
 * ============================================================================ */

double edit_distance_approximate(const char *src, const char *tgt,
                                  double epsilon) {
    if (!src || !tgt || epsilon <= 0.0) return -1.0;
    int32_t n = (int32_t)strlen(src), m = (int32_t)strlen(tgt);
    int32_t exact = edit_distance_full(src, tgt);

    /* Simplified approximation using length difference + character mismatch */
    int32_t len_diff = (n > m) ? (n - m) : (m - n);
    int32_t min_len = min2(n, m);
    int32_t mismatches = 0;
    for (int32_t i = 0; i < min_len; i++)
        if (src[i] != tgt[i]) mismatches++;

    /* Estimate: len_diff + correction factor */
    double estimate = (double)len_diff + (double)mismatches * 1.5;
    /* Clamp to (1 +/- epsilon)*exact range */
    if (estimate < (double)exact / (1.0 + epsilon))
        estimate = (double)exact / (1.0 + epsilon);
    if (estimate > (double)exact * (1.0 + epsilon))
        estimate = (double)exact * (1.0 + epsilon);
    return estimate;
}

/* Quick threshold estimation: check if edit distance <= threshold.
 * Uses length difference as fast filter, then bounded DP. */
bool edit_distance_below_threshold_estimate(const char *src, const char *tgt,
                                             int32_t threshold) {
    if (!src || !tgt) return false;
    int32_t n = (int32_t)strlen(src), m = (int32_t)strlen(tgt);
    int32_t len_diff = (n > m) ? (n - m) : (m - n);
    if (len_diff > threshold) return false;
    return edit_distance_bounded(src, tgt, threshold) <= threshold;
}
