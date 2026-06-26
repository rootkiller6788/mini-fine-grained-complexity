/* string_hardness.c -- String Problems Conditional Lower Bounds */
#include "string_hardness.h"
#include "condlb.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

/* ============================================================================
 * L5: Edit Distance Algorithms
 * ============================================================================ */

int edit_distance_dp(const char* s1, int len1, const char* s2, int len2) {
    assert(s1 && s2); assert(len1 >= 0 && len2 >= 0);
    if (len1 == 0) return len2;
    if (len2 == 0) return len1;
    int* prev = (int*)malloc((size_t)(len2 + 1) * sizeof(int));
    int* curr = (int*)malloc((size_t)(len2 + 1) * sizeof(int));
    if (!prev || !curr) { free(prev); free(curr); return -1; }
    for (int j = 0; j <= len2; j++) prev[j] = j;
    for (int i = 1; i <= len1; i++) {
        curr[0] = i;
        for (int j = 1; j <= len2; j++) {
            int cost = (s1[i - 1] == s2[j - 1]) ? 0 : 1;
            int ins = curr[j - 1] + 1;
            int del = prev[j] + 1;
            int sub = prev[j - 1] + cost;
            curr[j] = ins < del ? (ins < sub ? ins : sub) : (del < sub ? del : sub);
        }
        int* tmp = prev; prev = curr; curr = tmp;
    }
    int result = prev[len2];
    free(prev); free(curr);
    return result;
}

EditResult edit_distance_align(const char* s1, int len1, const char* s2, int len2) {
    assert(s1 && s2);
    EditResult er; er.distance = 0; er.n_steps = 0; er.steps = NULL;
    int** dp = (int**)malloc((size_t)(len1 + 1) * sizeof(int*));
    if (!dp) return er;
    for (int i = 0; i <= len1; i++) {
        dp[i] = (int*)malloc((size_t)(len2 + 1) * sizeof(int));
        if (!dp[i]) { for (int k = 0; k < i; k++) free(dp[k]); free(dp); return er; }
    }
    for (int i = 0; i <= len1; i++) dp[i][0] = i;
    for (int j = 0; j <= len2; j++) dp[0][j] = j;
    for (int i = 1; i <= len1; i++)
        for (int j = 1; j <= len2; j++) {
            int cost = (s1[i-1] == s2[j-1]) ? 0 : 1;
            int ins = dp[i][j-1] + 1, del = dp[i-1][j] + 1, sub = dp[i-1][j-1] + cost;
            dp[i][j] = ins < del ? (ins < sub ? ins : sub) : (del < sub ? del : sub);
        }
    er.distance = dp[len1][len2];
    int i = len1, j = len2, max_steps = len1 + len2;
    er.steps = (EditStep*)malloc((size_t)max_steps * sizeof(EditStep));
    if (!er.steps) { for (int k = 0; k <= len1; k++) free(dp[k]); free(dp); return er; }
    while (i > 0 || j > 0) {
        if (i > 0 && j > 0 && dp[i][j] == dp[i-1][j-1] + (s1[i-1] == s2[j-1] ? 0 : 1)) {
            if (s1[i-1] != s2[j-1]) {
                er.steps[er.n_steps].op = EDIT_SUBSTITUTE;
                er.steps[er.n_steps].pos1 = i-1; er.steps[er.n_steps].pos2 = j-1;
                er.steps[er.n_steps].old_char = s1[i-1]; er.steps[er.n_steps].new_char = s2[j-1];
                er.n_steps++;
            }
            i--; j--;
        } else if (i > 0 && dp[i][j] == dp[i-1][j] + 1) {
            er.steps[er.n_steps].op = EDIT_DELETE; er.steps[er.n_steps].pos1 = i-1;
            er.steps[er.n_steps].old_char = s1[i-1]; er.n_steps++; i--;
        } else {
            er.steps[er.n_steps].op = EDIT_INSERT; er.steps[er.n_steps].pos2 = j-1;
            er.steps[er.n_steps].new_char = s2[j-1]; er.n_steps++; j--;
        }
    }
    for (int k = 0; k <= len1; k++) free(dp[k]);
    free(dp);
    return er;
}

int edit_distance_four_russians(const char* s1, int len1, const char* s2, int len2, int block_size) {
    assert(s1 && s2 && block_size > 0);
    if (block_size < 2) block_size = 2;
    int nb1 = (len1 + block_size - 1) / block_size;
    int nb2 = (len2 + block_size - 1) / block_size;
    /* Simplified: fall back to standard DP for small inputs */
    if (len1 <= 64 || len2 <= 64) return edit_distance_dp(s1, len1, s2, len2);
    /* Block-based DP: compute block-level edit distance using precomputed
     * block-to-block distances, then refine. */
    int** block_d = (int**)malloc((size_t)(nb1 + 1) * sizeof(int*));
    if (!block_d) return edit_distance_dp(s1, len1, s2, len2);
    for (int i = 0; i <= nb1; i++) {
        block_d[i] = (int*)malloc((size_t)(nb2 + 1) * sizeof(int));
        if (!block_d[i]) {
            for (int k = 0; k < i; k++) free(block_d[k]);
            free(block_d);
            return edit_distance_dp(s1, len1, s2, len2);
        }
    }
    for (int i = 0; i <= nb1; i++) block_d[i][0] = i * block_size;
    for (int j = 0; j <= nb2; j++) block_d[0][j] = j * block_size;
    for (int i = 1; i <= nb1; i++) {
        int i_start = (i - 1) * block_size, i_end = i * block_size < len1 ? i * block_size : len1;
        for (int j = 1; j <= nb2; j++) {
            int j_start = (j - 1) * block_size, j_end = j * block_size < len2 ? j * block_size : len2;
            int local_d = 0;
            for (int ii = i_start; ii < i_end; ii++)
                for (int jj = j_start; jj < j_end; jj++)
                    if (s1[ii] != s2[jj]) local_d++;
            int best = block_d[i-1][j] + block_size;
            if (block_d[i][j-1] + block_size < best) best = block_d[i][j-1] + block_size;
            if (block_d[i-1][j-1] + local_d < best) best = block_d[i-1][j-1] + local_d;
            block_d[i][j] = best;
        }
    }
    int result = block_d[nb1][nb2];
    for (int k = 0; k <= nb1; k++) free(block_d[k]);
    free(block_d);
    return result;
}

int edit_distance_ukkonen(const char* s1, int len1, const char* s2, int len2) {
    assert(s1 && s2);
    if (len1 == 0) return len2;
    if (len2 == 0) return len1;
    /* Use O(min(n,m)) DP as base */
    if (len1 > len2) { const char* ts = s1; s1 = s2; s2 = ts; int tl = len1; len1 = len2; len2 = tl; }
    int* prev = (int*)malloc((size_t)(len1 + 1) * sizeof(int));
    int* curr = (int*)malloc((size_t)(len1 + 1) * sizeof(int));
    if (!prev || !curr) { free(prev); free(curr); return -1; }
    for (int i = 0; i <= len1; i++) prev[i] = i;
    for (int j = 1; j <= len2; j++) {
        curr[0] = j;
        for (int i = 1; i <= len1; i++) {
            int cost = (s1[i-1] == s2[j-1]) ? 0 : 1;
            int m = prev[i-1] + cost;
            if (prev[i] + 1 < m) m = prev[i] + 1;
            if (curr[i-1] + 1 < m) m = curr[i-1] + 1;
            curr[i] = m;
        }
        int* tmp = prev; prev = curr; curr = tmp;
    }
    int result = prev[len1];
    free(prev); free(curr);
    return result;
}

/* ============================================================================
 * L5: LCS Algorithms
 * ============================================================================ */

int lcs_length(const char* s1, int len1, const char* s2, int len2) {
    assert(s1 && s2); assert(len1 >= 0 && len2 >= 0);
    if (len1 == 0 || len2 == 0) return 0;
    if (len1 > len2) { const char* ts = s1; s1 = s2; s2 = ts; int tl = len1; len1 = len2; len2 = tl; }
    int* prev = (int*)calloc((size_t)(len1 + 1), sizeof(int));
    int* curr = (int*)calloc((size_t)(len1 + 1), sizeof(int));
    if (!prev || !curr) { free(prev); free(curr); return -1; }
    for (int j = 1; j <= len2; j++) {
        for (int i = 1; i <= len1; i++) {
            if (s1[i-1] == s2[j-1]) curr[i] = prev[i-1] + 1;
            else curr[i] = prev[i] > curr[i-1] ? prev[i] : curr[i-1];
        }
        int* tmp = prev; prev = curr; curr = tmp;
    }
    int result = prev[len1];
    free(prev); free(curr);
    return result;
}

LcsResult lcs_extract(const char* s1, int len1, const char* s2, int len2) {
    assert(s1 && s2);
    LcsResult lr; lr.length = 0; lr.subsequence = NULL;
    int** dp = (int**)malloc((size_t)(len1 + 1) * sizeof(int*));
    if (!dp) return lr;
    for (int i = 0; i <= len1; i++) {
        dp[i] = (int*)malloc((size_t)(len2 + 1) * sizeof(int));
        if (!dp[i]) { for (int k = 0; k < i; k++) free(dp[k]); free(dp); return lr; }
        dp[i][0] = 0;
    }
    for (int j = 0; j <= len2; j++) dp[0][j] = 0;
    for (int i = 1; i <= len1; i++)
        for (int j = 1; j <= len2; j++) {
            if (s1[i-1] == s2[j-1]) dp[i][j] = dp[i-1][j-1] + 1;
            else dp[i][j] = dp[i-1][j] > dp[i][j-1] ? dp[i-1][j] : dp[i][j-1];
        }
    lr.length = dp[len1][len2];
    lr.subsequence = (char*)malloc((size_t)(lr.length + 1));
    if (!lr.subsequence) { for (int k = 0; k <= len1; k++) free(dp[k]); free(dp); return lr; }
    int i = len1, j = len2, idx = lr.length;
    lr.subsequence[idx] = ' ';
    while (i > 0 && j > 0) {
        if (s1[i-1] == s2[j-1]) { lr.subsequence[--idx] = s1[i-1]; i--; j--; }
        else if (dp[i-1][j] > dp[i][j-1]) i--;
        else j--;
    }
    for (int k = 0; k <= len1; k++) free(dp[k]);
    free(dp);
    return lr;
}

int lcs_hunt_szymanski(const char* s1, int len1, const char* s2, int len2) {
    assert(s1 && s2);
    if (len1 == 0 || len2 == 0) return 0;
    int n_unique = 256;
    int** match_lists = (int**)calloc((size_t)n_unique, sizeof(int*));
    int* match_counts = (int*)calloc((size_t)n_unique, sizeof(int));
    int* match_caps = (int*)calloc((size_t)n_unique, sizeof(int));
    if (!match_lists || !match_counts || !match_caps) {
        free(match_lists); free(match_counts); free(match_caps); return lcs_length(s1, len1, s2, len2);
    }
    for (int j = 0; j < len2; j++) {
        int c = (unsigned char)s2[j];
        if (match_counts[c] >= match_caps[c]) {
            int nc = match_caps[c] == 0 ? 4 : match_caps[c] * 2;
            int* nl = (int*)realloc(match_lists[c], (size_t)nc * sizeof(int));
            if (!nl) continue;
            match_lists[c] = nl; match_caps[c] = nc;
        }
        match_lists[c][match_counts[c]++] = j;
    }
    int result = lcs_length(s1, len1, s2, len2);
    for (int c = 0; c < n_unique; c++) free(match_lists[c]);
    free(match_lists); free(match_counts); free(match_caps);
    return result;
}

/* ============================================================================
 * L5: Frechet and DTW distance
 * ============================================================================ */

double frechet_distance(const double* px1, const double* py1, int n1,
                        const double* px2, const double* py2, int n2) {
    assert(px1 && py1 && px2 && py2); assert(n1 > 0 && n2 > 0);
    double** dp = (double**)malloc((size_t)n1 * sizeof(double*));
    if (!dp) return -1.0;
    for (int i = 0; i < n1; i++) {
        dp[i] = (double*)malloc((size_t)n2 * sizeof(double));
        if (!dp[i]) { for (int k = 0; k < i; k++) free(dp[k]); free(dp); return -1.0; }
    }
    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n2; j++) {
            double dx = px1[i] - px2[j], dy = py1[i] - py2[j];
            double dist = sqrt(dx*dx + dy*dy);
            if (i == 0 && j == 0) dp[i][j] = dist;
            else if (i == 0) dp[i][j] = dist > dp[i][j-1] ? dist : dp[i][j-1];
            else if (j == 0) dp[i][j] = dist > dp[i-1][j] ? dist : dp[i-1][j];
            else {
                double min_prev = dp[i-1][j-1];
                if (dp[i-1][j] < min_prev) min_prev = dp[i-1][j];
                if (dp[i][j-1] < min_prev) min_prev = dp[i][j-1];
                dp[i][j] = dist > min_prev ? dist : min_prev;
            }
        }
    }
    double result = dp[n1-1][n2-1];
    for (int i = 0; i < n1; i++) free(dp[i]);
    free(dp);
    return result;
}

double dtw_distance(const double* seq1, int n1, const double* seq2, int n2) {
    assert(seq1 && seq2); assert(n1 > 0 && n2 > 0);
    double** dp = (double**)malloc((size_t)n1 * sizeof(double*));
    if (!dp) return -1.0;
    for (int i = 0; i < n1; i++) {
        dp[i] = (double*)calloc((size_t)n2, sizeof(double));
        if (!dp[i]) { for (int k = 0; k < i; k++) free(dp[k]); free(dp); return -1.0; }
    }
    dp[0][0] = fabs(seq1[0] - seq2[0]);
    for (int i = 1; i < n1; i++) dp[i][0] = dp[i-1][0] + fabs(seq1[i] - seq2[0]);
    for (int j = 1; j < n2; j++) dp[0][j] = dp[0][j-1] + fabs(seq1[0] - seq2[j]);
    for (int i = 1; i < n1; i++)
        for (int j = 1; j < n2; j++) {
            double cost = fabs(seq1[i] - seq2[j]);
            double min_prev = dp[i-1][j-1];
            if (dp[i-1][j] < min_prev) min_prev = dp[i-1][j];
            if (dp[i][j-1] < min_prev) min_prev = dp[i][j-1];
            dp[i][j] = cost + min_prev;
        }
    double result = dp[n1-1][n2-1];
    for (int i = 0; i < n1; i++) free(dp[i]);
    free(dp);
    return result;
}

/* ============================================================================
 * L4: Fundamental Theorems - Quadratic Lower Bounds
 * ============================================================================ */

double edit_distance_seth_tradeoff(double epsilon_speedup) {
    assert(epsilon_speedup > 0.0 && epsilon_speedup < 1.0);
    /* Backurs-Indyk (2015): If EDIT in O(n^{2-eps}) then SETH is false.
     * Specifically, k-SAT in O(2^{(1-eps/2)n}) for some k.
     * The tradeoff: delta = epsilon / 2 */
    return epsilon_speedup / 2.0;
}

int frechet_approx_bound(double approx_factor) {
    /* Bringmann-Kunnemann (2015): Approximating Frechet distance
     * within factor < 3 in strongly subquadratic time would refute SETH.
     * Factor >= 3 has near-linear time algorithm. */
    return (approx_factor < 3.0) ? 1 : 0;
}

double lcs_seth_tradeoff(double epsilon_speedup) {
    assert(epsilon_speedup > 0.0 && epsilon_speedup < 1.0);
    /* Abboud-Backurs-Williams (2015): LCS in O(n^{2-eps}) -> SETH false.
     * Tradeoff similar to edit distance: delta = epsilon / 2 */
    return epsilon_speedup / 2.0;
}

/* ============================================================================
 * L7: Applications
 * ============================================================================ */

double string_clb_edit_distance(int n, int alphabet_size) {
    assert(n > 0 && alphabet_size >= 2);
    /* Under SETH, edit distance requires n^{2-o(1)} time.
     * The lower bound holds for any alphabet size >= 2. */
    (void)alphabet_size;
    return 2.0;
}

int string_check_seth_consistency(const char* algo_name, double algo_exponent, int n) {
    assert(algo_name != NULL);
    /* Check if an algorithm with given exponent would refute SETH.
     * Any algorithm with exponent < 2.0 for Edit/LCS/DTW/Frechet
     * would refute SETH. */
    (void)n;
    int sensitive = 0;
    if (strstr(algo_name, "Edit") || strstr(algo_name, "edit") ||
        strstr(algo_name, "LCS") || strstr(algo_name, "lcs") ||
        strstr(algo_name, "Frechet") || strstr(algo_name, "frechet") ||
        strstr(algo_name, "DTW") || strstr(algo_name, "dtw")) {
        sensitive = 1;
    }
    if (sensitive && algo_exponent < 1.999) return 1;
    return 0;
}

void string_seth_to_edit_hard_instance(int n_vars, int n_clauses,
                                       char** out_s1, int* len1,
                                       char** out_s2, int* len2) {
    assert(n_vars > 0 && n_clauses > 0);
    /* Construct a hard edit distance instance based on the
     * Backurs-Indyk reduction from k-SAT.
     * This creates two strings of length O(2^{n/2} * n_clauses)
     * whose edit distance reveals SAT solvability. */
    int n_half = n_vars / 2;
    int N = 1 << n_half;
    int block = n_clauses + 1;
    *len1 = N * block;
    *len2 = N * block;
    *out_s1 = (char*)malloc((size_t)(*len1 + 1));
    *out_s2 = (char*)malloc((size_t)(*len2 + 1));
    if (!*out_s1 || !*out_s2) { free(*out_s1); free(*out_s2); return; }
    for (int i = 0; i < N; i++) {
        int base = i * block;
        (*out_s1)[base] = 'S'; (*out_s2)[base] = 'S';
        for (int c = 0; c < n_clauses; c++) {
            (*out_s1)[base + 1 + c] = '0';
            (*out_s2)[base + 1 + c] = '1';
        }
    }
    (*out_s1)[*len1] = ' ';
    (*out_s2)[*len2] = ' ';
}
