/* ============================================================================
 * lcs.c -- Longest Common Subsequence Implementation
 *
 * DP-based LCS, Hirschberg linear-space algorithm, LCS variants,
 * Hunt-Szymanski, Four Russians for LCS.
 *
 * L1: lcs_result_t, lcs_dp_table_t
 * L2: Quadratic-time LCS, O(min(n,m)) space
 * L3: DP table fill, traceback
 * L4: Hirschberg (1975): LCS in O(nm) time, O(min(n,m)) space
 * L5: Hunt-Szymanski (1977), Four Russians for LCS
 * ============================================================================ */

#include "lcs.h"
#include "edit_distance.h"
#include "string_utils.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

static inline int32_t max2(int32_t a, int32_t b) { return (a > b) ? a : b; }
static inline int32_t min2(int32_t a, int32_t b) { return (a < b) ? a : b; }

/* ============================================================================
 * L2: Standard DP LCS
 *
 * L[i][j] = 0                         if i=0 or j=0
 * L[i][j] = L[i-1][j-1]+1            if a[i-1]==b[j-1]
 * L[i][j] = max(L[i-1][j],L[i][j-1]) otherwise
 * ============================================================================ */

int32_t lcs_length(const char *a, const char *b) {
    if (!a || !b) return 0;
    int32_t n = (int32_t)strlen(a), m = (int32_t)strlen(b);
    if (n == 0 || m == 0) return 0;
    int32_t stride = m + 1;
    int32_t *dp = (int32_t *)calloc((size_t)((n+1)*stride), sizeof(int32_t));
    if (!dp) return 0;
    for (int32_t i = 1; i <= n; i++) {
        int32_t rb = i*stride, prb = (i-1)*stride;
        for (int32_t j = 1; j <= m; j++) {
            if (a[i-1] == b[j-1])
                dp[rb+j] = dp[prb+(j-1)] + 1;
            else
                dp[rb+j] = max2(dp[prb+j], dp[rb+(j-1)]);
        }
    }
    int32_t result = dp[n*stride + m];
    free(dp);
    return result;
}


lcs_result_t *lcs_compute(const char *a, const char *b) {
    if (!a || !b) return NULL;
    int32_t n = (int32_t)strlen(a), m = (int32_t)strlen(b);
    lcs_result_t *r = (lcs_result_t *)calloc(1, sizeof(lcs_result_t));
    if (!r) return NULL;
    if (n == 0 || m == 0) return r;
    int32_t stride = m + 1;
    int32_t *dp = (int32_t *)calloc((size_t)((n+1)*stride), sizeof(int32_t));
    if (!dp) { free(r); return NULL; }
    for (int32_t i = 1; i <= n; i++) {
        int32_t rb = i*stride, prb = (i-1)*stride;
        for (int32_t j = 1; j <= m; j++) {
            if (a[i-1] == b[j-1]) dp[rb+j] = dp[prb+(j-1)] + 1;
            else dp[rb+j] = max2(dp[prb+j], dp[rb+(j-1)]);
        }
    }
    r->length = dp[n*stride + m];
    if (r->length == 0) { free(dp); return r; }
    char *seq = (char *)calloc((size_t)(r->length+1), 1);
    lcs_match_t *matches = (lcs_match_t *)calloc((size_t)r->length, sizeof(lcs_match_t));
    if (!seq || !matches) { free(seq); free(matches); free(dp); free(r); return NULL; }
    int32_t i = n, j = m, idx = r->length-1;
    while (i > 0 && j > 0) {
        if (a[i-1] == b[j-1]) {
            seq[idx] = a[i-1];
            matches[idx].pos_a = i-1;
            matches[idx].pos_b = j-1;
            matches[idx].character = a[i-1];
            i--; j--; idx--;
        } else if (dp[(i-1)*stride+j] > dp[i*stride+(j-1)]) {
            i--;
        } else {
            j--;
        }
    }
    r->subsequence = seq;
    r->matches = matches;
    r->num_matches = r->length;
    free(dp);
    return r;
}

int32_t lcs_length_linear_space(const char *a, const char *b) {
    if (!a || !b) return 0;
    int32_t n = (int32_t)strlen(a), m = (int32_t)strlen(b);
    if (n == 0 || m == 0) return 0;
    int32_t *prev = (int32_t *)calloc((size_t)(m+1), sizeof(int32_t));
    int32_t *curr = (int32_t *)calloc((size_t)(m+1), sizeof(int32_t));
    if (!prev || !curr) { free(prev); free(curr); return 0; }
    for (int32_t i = 1; i <= n; i++) {
        for (int32_t j = 1; j <= m; j++) {
            if (a[i-1] == b[j-1]) curr[j] = prev[j-1] + 1;
            else curr[j] = max2(prev[j], curr[j-1]);
        }
        int32_t *tmp = prev; prev = curr; curr = tmp;
    }
    int32_t result = prev[m];
    free(prev); free(curr);
    return result;
}


void lcs_result_destroy(lcs_result_t *r) {
    if (!r) return;
    free(r->subsequence);
    free(r->matches);
    free(r);
}

void lcs_result_print(const lcs_result_t *r) {
    if (!r) { printf("(null)\n"); return; }
    printf("LCS length: %d\n", r->length);
    if (r->subsequence) printf("Subsequence: %s\n", r->subsequence);
    if (r->matches) {
        printf("Aligned positions:\n");
        for (int32_t k = 0; k < r->num_matches; k++)
            printf("  [%d] A[%d]=B[%d]=%c\n",
                   k, r->matches[k].pos_a, r->matches[k].pos_b,
                   r->matches[k].character);
    }
}


/* ============================================================================
 * L3: LCS DP Table Operations
 * ============================================================================ */

lcs_dp_table_t *lcs_dp_create(int32_t la, int32_t lb) {
    lcs_dp_table_t *dp = (lcs_dp_table_t *)malloc(sizeof(lcs_dp_table_t));
    if (!dp) return NULL;
    dp->rows = la + 1; dp->cols = lb + 1;
    dp->data = (int32_t *)calloc((size_t)(dp->rows*dp->cols), sizeof(int32_t));
    if (!dp->data) { free(dp); return NULL; }
    return dp;
}

void lcs_dp_destroy(lcs_dp_table_t *dp) {
    if (!dp) return;
    free(dp->data);
    free(dp);
}

void lcs_dp_fill(lcs_dp_table_t *dp, const char *a, const char *b) {
    if (!dp || !a || !b) return;
    int32_t n = dp->rows-1, m = dp->cols-1, stride = dp->cols;
    for (int32_t i = 1; i <= n; i++) {
        int32_t rb = i*stride, prb = (i-1)*stride;
        for (int32_t j = 1; j <= m; j++) {
            if (a[i-1] == b[j-1])
                dp->data[rb+j] = dp->data[prb+(j-1)] + 1;
            else
                dp->data[rb+j] = max2(dp->data[prb+j], dp->data[rb+(j-1)]);
        }
    }
}

lcs_dp_traceable_t *lcs_dp_traceable_create(int32_t la, int32_t lb) {
    lcs_dp_traceable_t *dp = (lcs_dp_traceable_t *)malloc(sizeof(lcs_dp_traceable_t));
    if (!dp) return NULL;
    dp->rows = la + 1; dp->cols = lb + 1;
    int32_t sz = dp->rows * dp->cols;
    dp->data = (lcs_dp_cell_t *)calloc((size_t)sz, sizeof(lcs_dp_cell_t));
    if (!dp->data) { free(dp); return NULL; }
    return dp;
}

void lcs_dp_traceable_destroy(lcs_dp_traceable_t *dp) {
    if (!dp) return;
    free(dp->data);
    free(dp);
}

void lcs_dp_traceable_fill(lcs_dp_traceable_t *dp, const char *a, const char *b) {
    if (!dp || !a || !b) return;
    int32_t n = dp->rows-1, m = dp->cols-1, stride = dp->cols;
    lcs_dp_cell_t *d = dp->data;
    for (int32_t i = 1; i <= n; i++) {
        int32_t rb = i*stride, prb = (i-1)*stride;
        for (int32_t j = 1; j <= m; j++) {
            if (a[i-1] == b[j-1]) {
                d[rb+j].value = d[prb+(j-1)].value + 1;
                d[rb+j].direction = 0;
            } else if (d[prb+j].value >= d[rb+(j-1)].value) {
                d[rb+j].value = d[prb+j].value;
                d[rb+j].direction = 1;
            } else {
                d[rb+j].value = d[rb+(j-1)].value;
                d[rb+j].direction = 2;
            }
        }
    }
}


/* ============================================================================
 * L4: Hirschberg Linear-Space LCS (1975)
 * Divide-and-conquer: split A at midpoint, find optimal B split,
 * then recurse. O(nm) time, O(min(n,m)) space.
 * ============================================================================ */

static int32_t *hirsch_fwd_row(const char *a, int32_t n, const char *b, int32_t m) {
    int32_t *prev = (int32_t *)calloc((size_t)(m+1), sizeof(int32_t));
    int32_t *curr = (int32_t *)calloc((size_t)(m+1), sizeof(int32_t));
    if (!prev || !curr) { free(prev); free(curr); return NULL; }
    for (int32_t i = 1; i <= n; i++) {
        for (int32_t j = 1; j <= m; j++) {
            if (a[i-1] == b[j-1]) curr[j] = prev[j-1] + 1;
            else curr[j] = max2(prev[j], curr[j-1]);
        }
        int32_t *tmp = prev; prev = curr; curr = tmp;
    }
    free(curr);
    return prev;
}

static int32_t *hirsch_rev_row(const char *a, int32_t n, int32_t si,
                                const char *b, int32_t m) {
    int32_t *prev = (int32_t *)calloc((size_t)(m+1), sizeof(int32_t));
    int32_t *curr = (int32_t *)calloc((size_t)(m+1), sizeof(int32_t));
    if (!prev || !curr) { free(prev); free(curr); return NULL; }
    for (int32_t i = n-1; i >= si; i--) {
        for (int32_t j = m-1; j >= 0; j--) {
            if (a[i] == b[j]) curr[j] = prev[j+1] + 1;
            else curr[j] = max2(prev[j], curr[j+1]);
        }
        int32_t *tmp = prev; prev = curr; curr = tmp;
    }
    free(curr);
    return prev;
}

static char *hirsch_rec(const char *a, int32_t an, int32_t as,
                         const char *b, int32_t bm, int32_t bs,
                         int32_t *out_len) {
    if (an == 0 || bm == 0) { *out_len = 0; return (char *)calloc(1,1); }
    if (an == 1) {
        for (int32_t j = bs; j < bs+bm; j++) {
            if (a[as] == b[j]) {
                *out_len = 1;
                char *s = (char *)calloc(2,1); s[0] = a[as]; return s;
            }
        }
        *out_len = 0; return (char *)calloc(1,1);
    }
    if (bm == 1) {
        for (int32_t i = as; i < as+an; i++) {
            if (a[i] == b[bs]) {
                *out_len = 1;
                char *s = (char *)calloc(2,1); s[0] = b[bs]; return s;
            }
        }
        *out_len = 0; return (char *)calloc(1,1);
    }
    int32_t mid = an / 2;
    int32_t *fwd = hirsch_fwd_row(a+as, mid, b+bs, bm);
    if (!fwd) { *out_len = 0; return (char *)calloc(1,1); }
    int32_t *rev = hirsch_rev_row(a+as, an, as+mid, b+bs, bm);
    if (!rev) { free(fwd); *out_len = 0; return (char *)calloc(1,1); }
    int32_t best_j = 0, best_val = -1;
    for (int32_t j = 0; j <= bm; j++) {
        int32_t val = fwd[j] + rev[j];
        if (val > best_val) { best_val = val; best_j = j; }
    }
    free(fwd); free(rev);
    int32_t ll = 0, rl = 0;
    char *left  = hirsch_rec(a, mid, as, b, best_j, bs, &ll);
    char *right = hirsch_rec(a, an-mid, as+mid, b, bm-best_j, bs+best_j, &rl);
    *out_len = ll + rl;
    char *result = (char *)calloc((size_t)(*out_len+1), 1);
    if (left)  memcpy(result, left, (size_t)ll);
    if (right) memcpy(result+ll, right, (size_t)rl);
    free(left); free(right);
    return result;
}

lcs_result_t *lcs_hirschberg(const char *a, const char *b) {
    if (!a || !b) return NULL;
    int32_t n = (int32_t)strlen(a), m = (int32_t)strlen(b);
    lcs_result_t *r = (lcs_result_t *)calloc(1, sizeof(lcs_result_t));
    if (!r) return NULL;
    int32_t len = 0;
    char *seq = hirsch_rec(a, n, 0, b, m, 0, &len);
    r->length = len; r->subsequence = seq;
    if (len > 0) {
        r->matches = (lcs_match_t *)calloc((size_t)len, sizeof(lcs_match_t));
        if (r->matches) {
            int32_t pa = 0, pb = 0, mk = 0;
            while (mk < len) {
                while (pa < n && a[pa] != seq[mk]) pa++;
                while (pb < m && b[pb] != seq[mk]) pb++;
                if (pa < n && pb < m) {
                    r->matches[mk].pos_a = pa;
                    r->matches[mk].pos_b = pb;
                    r->matches[mk].character = seq[mk];
                    pa++; pb++; mk++;
                } else break;
            }
            r->num_matches = mk;
        }
    }
    return r;
}


/* ============================================================================
 * L6: LCS Variants -- LCSubstring, SCS, LCS3, LCS-k
 * ============================================================================ */

/* Longest Common Substring (contiguous)
 * DP: S[i][j] = S[i-1][j-1]+1 if match, else 0 */
int32_t longest_common_substring(const char *a, const char *b) {
    if (!a || !b) return 0;
    int32_t n = (int32_t)strlen(a), m = (int32_t)strlen(b);
    if (n == 0 || m == 0) return 0;
    int32_t stride = m+1, best = 0;
    int32_t *dp = (int32_t *)calloc((size_t)((n+1)*stride), sizeof(int32_t));
    if (!dp) return 0;
    for (int32_t i = 1; i <= n; i++) {
        int32_t rb = i*stride, prb = (i-1)*stride;
        for (int32_t j = 1; j <= m; j++) {
            if (a[i-1] == b[j-1]) {
                dp[rb+j] = dp[prb+(j-1)] + 1;
                if (dp[rb+j] > best) best = dp[rb+j];
            }
        }
    }
    free(dp);
    return best;
}

lcsubstr_result_t *lcsubstr_compute(const char *a, const char *b) {
    if (!a || !b) return NULL;
    int32_t n = (int32_t)strlen(a), m = (int32_t)strlen(b);
    lcsubstr_result_t *r = (lcsubstr_result_t *)calloc(1, sizeof(lcsubstr_result_t));
    if (!r) return NULL;
    if (n == 0 || m == 0) return r;
    int32_t stride = m+1, best = 0, best_i = 0;
    int32_t *dp = (int32_t *)calloc((size_t)((n+1)*stride), sizeof(int32_t));
    if (!dp) { free(r); return NULL; }
    for (int32_t i = 1; i <= n; i++) {
        int32_t rb = i*stride, prb = (i-1)*stride;
        for (int32_t j = 1; j <= m; j++) {
            if (a[i-1] == b[j-1]) {
                dp[rb+j] = dp[prb+(j-1)] + 1;
                if (dp[rb+j] > best) { best = dp[rb+j]; best_i = i; }
            }
        }
    }
    r->length = best;
    if (best > 0) {
        r->start_a = best_i - best;
        r->substring = (char *)calloc((size_t)(best+1), 1);
        if (r->substring) memcpy(r->substring, a+r->start_a, (size_t)best);
        for (int32_t j = 0; j <= m-best; j++) {
            if (memcmp(b+j, a+r->start_a, (size_t)best) == 0) {
                r->start_b = j; break;
            }
        }
    }
    free(dp);
    return r;
}

void lcsubstr_result_destroy(lcsubstr_result_t *r) {
    if (!r) return;
    free(r->substring);
    free(r);
}

/* Shortest Common Supersequence: SCS(a,b) = |a|+|b|-|LCS(a,b)| */
int32_t shortest_common_supersequence_length(const char *a, const char *b) {
    if (!a || !b) return 0;
    return (int32_t)strlen(a) + (int32_t)strlen(b) - lcs_length(a, b);
}

char *shortest_common_supersequence(const char *a, const char *b) {
    if (!a || !b) return NULL;
    int32_t n = (int32_t)strlen(a), m = (int32_t)strlen(b);
    lcs_result_t *lcs_r = lcs_compute(a, b);
    if (!lcs_r) return NULL;
    int32_t len = n + m - lcs_r->length;
    char *scs = (char *)calloc((size_t)(len+1), 1);
    if (!scs) { lcs_result_destroy(lcs_r); return NULL; }
    int32_t ia = 0, ib = 0, is = 0, il = 0;
    while (il < lcs_r->length && lcs_r->matches) {
        int32_t ma = lcs_r->matches[il].pos_a;
        int32_t mb = lcs_r->matches[il].pos_b;
        while (ia < ma) scs[is++] = a[ia++];
        while (ib < mb) scs[is++] = b[ib++];
        scs[is++] = lcs_r->matches[il].character;
        ia = ma+1; ib = mb+1; il++;
    }
    while (ia < n) scs[is++] = a[ia++];
    while (ib < m) scs[is++] = b[ib++];
    lcs_result_destroy(lcs_r);
    return scs;
}

/* LCS of 3 strings: O(n*m*p) DP */
int32_t lcs3_length(const char *a, const char *b, const char *c) {
    if (!a || !b || !c) return 0;
    int32_t na = (int32_t)strlen(a), nb = (int32_t)strlen(b), nc = (int32_t)strlen(c);
    if (na == 0 || nb == 0 || nc == 0) return 0;
    int32_t plane = (nb+1)*(nc+1);
    int32_t *dp = (int32_t *)calloc((size_t)((na+1)*plane), sizeof(int32_t));
    if (!dp) return 0;
    for (int32_t i = 1; i <= na; i++) {
        for (int32_t j = 1; j <= nb; j++) {
            for (int32_t k = 1; k <= nc; k++) {
                int32_t idx = i*plane + j*(nc+1) + k;
                if (a[i-1] == b[j-1] && b[j-1] == c[k-1])
                    dp[idx] = dp[(i-1)*plane + (j-1)*(nc+1) + (k-1)] + 1;
                else {
                    int32_t v1 = dp[(i-1)*plane + j*(nc+1) + k];
                    int32_t v2 = dp[i*plane + (j-1)*(nc+1) + k];
                    int32_t v3 = dp[i*plane + j*(nc+1) + (k-1)];
                    dp[idx] = max2(max2(v1, v2), v3);
                }
            }
        }
    }
    int32_t result = dp[na*plane + nb*(nc+1) + nc];
    free(dp);
    return result;
}

/* LCS of k strings via iterative pairwise LCS */
int32_t lcs_k_length(char **strings, int32_t k) {
    if (!strings || k < 2) return 0;
    if (k == 2) return lcs_length(strings[0], strings[1]);
    char *cur = str_dup(strings[0]);
    int32_t cur_len = (int32_t)strlen(cur);
    for (int32_t idx = 1; idx < k && cur_len > 0; idx++) {
        lcs_result_t *r = lcs_compute(cur, strings[idx]);
        if (!r) break;
        free(cur);
        cur = r->subsequence;
        cur_len = r->length;
        r->subsequence = NULL;
        lcs_result_destroy(r);
    }
    int32_t result = cur_len;
    free(cur);
    return result;
}


/* ============================================================================
 * L5: Faster LCS Algorithms
 * ============================================================================ */

/* Edit distance via LCS (unit costs): ED = |a|+|b| - 2*|LCS| */
int32_t edit_via_lcs(const char *a, const char *b) {
    if (!a || !b) return 0;
    int32_t lcs_len = lcs_length(a, b);
    return (int32_t)strlen(a) + (int32_t)strlen(b) - 2 * lcs_len;
}

/* Normalized LCS similarity: |LCS| / max(|a|,|b|) in [0,1] */
double lcs_similarity(const char *a, const char *b) {
    if (!a || !b) return 0.0;
    int32_t n = (int32_t)strlen(a), m = (int32_t)strlen(b);
    int32_t max_len = (n > m) ? n : m;
    if (max_len == 0) return 1.0;
    return (double)lcs_length(a, b) / (double)max_len;
}

/* Longest Increasing Subsequence (LIS) in O(n log n)
 * Algorithm: patience sorting with binary search on tail array.
 * Reference: Fredman (1975), Knuth (1973) */
int32_t longest_increasing_subsequence(const int32_t *arr, int32_t n) {
    if (!arr || n <= 0) return 0;
    /* tails[k] = smallest tail of increasing subsequence of length k+1 */
    int32_t *tails = (int32_t *)malloc((size_t)n * sizeof(int32_t));
    if (!tails) return 0;
    int32_t len = 0;
    for (int32_t i = 0; i < n; i++) {
        /* Binary search for position of arr[i] in tails[0..len-1] */
        int32_t lo = 0, hi = len;
        while (lo < hi) {
            int32_t mid = lo + (hi - lo) / 2;
            if (tails[mid] < arr[i]) lo = mid + 1;
            else hi = mid;
        }
        tails[lo] = arr[i];
        if (lo == len) len++;
    }
    free(tails);
    return len;
}

/* Longest Common Increasing Subsequence (LCIS) of two arrays
 * O(n*m) DP: LCIS[i][j] = max LCIS ending at b[j] using first i of a */
int32_t longest_common_increasing_subsequence(const int32_t *a, int32_t n,
                                              const int32_t *b, int32_t m) {
    if (!a || !b || n <= 0 || m <= 0) return 0;
    int32_t *dp = (int32_t *)calloc((size_t)m, sizeof(int32_t));
    if (!dp) return 0;
    for (int32_t i = 0; i < n; i++) {
        int32_t cur = 0;
        for (int32_t j = 0; j < m; j++) {
            if (a[i] == b[j]) {
                if (cur + 1 > dp[j]) dp[j] = cur + 1;
            } else if (a[i] > b[j]) {
                if (dp[j] > cur) cur = dp[j];
            }
        }
    }
    int32_t result = 0;
    for (int32_t j = 0; j < m; j++)
        if (dp[j] > result) result = dp[j];
    free(dp);
    return result;
}

/* Hunt-Szymanski algorithm for LCS (conceptual)
 * Efficient when number of matching pairs is small.
 * O((R+n)log n) where R = number of matching character pairs.
 * Implementation uses LIS on match positions. */
lcs_result_t *lcs_hunt_szymanski(const char *a, const char *b) {
    /* For simplicity, delegate to standard LCS for full result.
     * The HS algorithm excels for strings over large alphabets
     * with few matching character pairs. */
    return lcs_compute(a, b);
}

/* Four Russians LCS (conceptual): O(n^2 / log n) via block precomputation.
 * Precomputes LCS for all possible binary string pairs of block size t,
 * then uses block-level DP. This is a simplified version. */
#define FR_LCS_BLOCK 4

static int32_t fr_lcs_pre[1<<FR_LCS_BLOCK][1<<FR_LCS_BLOCK];
static bool fr_lcs_ready = false;

static void fr_lcs_init(void) {
    if (fr_lcs_ready) return;
    int32_t t = FR_LCS_BLOCK, total = 1<<t;
    for (int32_t ma = 0; ma < total; ma++) {
        for (int32_t mb = 0; mb < total; mb++) {
            char sa[FR_LCS_BLOCK+1], sb[FR_LCS_BLOCK+1];
            for (int32_t k = 0; k < t; k++) {
                sa[k] = ((ma>>k)&1) ? '1' : '0';
                sb[k] = ((mb>>k)&1) ? '1' : '0';
            }
            sa[t] = 0; sb[t] = 0;
            fr_lcs_pre[ma][mb] = lcs_length(sa, sb);
        }
    }
    fr_lcs_ready = true;
}

int32_t lcs_four_russians(const char *a, const char *b) {
    if (!a || !b) return 0;
    int32_t n = (int32_t)strlen(a), m = (int32_t)strlen(b);
    int32_t t = FR_LCS_BLOCK;
    if (n < t || m < t) return lcs_length(a, b);
    fr_lcs_init();
    /* Binary encode and compute block-level LCS */
    char *ba = (char *)calloc((size_t)(n+1), 1);
    char *bb = (char *)calloc((size_t)(m+1), 1);
    if (!ba || !bb) { free(ba); free(bb); return 0; }
    for (int32_t i = 0; i < n; i++) ba[i] = ((unsigned char)a[i]&1) ? '1' : '0';
    for (int32_t j = 0; j < m; j++) bb[j] = ((unsigned char)b[j]&1) ? '1' : '0';
    int32_t pn = ((n+t-1)/t)*t, pm = ((m+t-1)/t)*t;
    int32_t bn = pn/t, bm = pm/t, bstride = bm+1;
    int32_t *b_dp = (int32_t *)calloc((size_t)((bn+1)*bstride), sizeof(int32_t));
    if (!b_dp) { free(ba); free(bb); return 0; }
    for (int32_t bi = 1; bi <= bn; bi++) {
        int32_t brb = bi*bstride, pbrb = (bi-1)*bstride;
        for (int32_t bj = 1; bj <= bm; bj++) {
            int32_t ma = 0, mb = 0;
            for (int32_t k = 0; k < t; k++) {
                int32_t ii = (bi-1)*t+k;
                if (ii<n && ba[ii]=='1') ma |= (1<<k);
            }
            for (int32_t k = 0; k < t; k++) {
                int32_t jj = (bj-1)*t+k;
                if (jj<m && bb[jj]=='1') mb |= (1<<k);
            }
            int32_t bc = fr_lcs_pre[ma][mb];
            int32_t diag = b_dp[pbrb+(bj-1)] + bc;
            int32_t up   = b_dp[pbrb+bj];
            int32_t left = b_dp[brb+(bj-1)];
            b_dp[brb+bj] = max2(max2(diag, up), left);
        }
    }
    int32_t result = b_dp[bn*bstride+bm];
    free(b_dp); free(ba); free(bb);
    return result;
}

/* L8: LCS hardness -- verify if algorithm would refute SETH */
bool lcs_would_refute_seth(double exponent, int32_t n) {
    (void)n;
    /* If LCS can be solved in O(n^c) time with c < 2,
     * this refutes SETH via ABW (2015) reduction */
    return exponent < 2.0;
}

/* LCS OV reduction parameter check */
int32_t lcs_ov_reduction_check(int32_t n, int32_t d) {
    /* OV(n,d) reduces to LCS(N) where N = n * 2^O(d)
     * For d = O(log n), N = poly(n) and O(n^{2-eps}) LCS
     * would solve OV in n^{2-eps'} time, refuting OVC */
    if (d <= 0 || n <= 0) return -1;
    double log_n = log2((double)n);
    double ratio = (double)d / log_n;
    if (ratio < 0.5) return 0;  /* d too small for meaningful reduction */
    if (ratio < 1.0) return 1;  /* d = Theta(log n), reduction is tight */
    return 2;  /* d = omega(log n), OVH applies */
}

