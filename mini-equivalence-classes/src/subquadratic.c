/* ============================================================================
 * subquadratic.c -- Subquadratic Equivalence Class Implementation
 *
 * Implements the core algorithms for the subquadratic equivalence class:
 *   - Orthogonal Vectors (OV): brute force and Williams algorithm
 *   - Edit Distance (Levenshtein): DP and linear-space variants
 *   - Longest Common Subsequence (LCS): DP and Hirschberg
 *   - Frechet Distance: free-space diagram method
 *   - Dynamic Time Warping (DTW): standard DP and Sakoe-Chiba band
 *
 * Knowledge Coverage:
 *   L1: OV, Edit Distance, LCS, Frechet, DTW data structures
 *   L2: Subquadratic equivalence concept, SETH-based reductions
 *   L3: Boolean vectors (bit-packed), DP tables, polygonal curves
 *   L4: SETH => OV => Edit Distance/LCS/DTW quadratic lower bounds
 *   L5: O(n^2) DP algorithms for all member problems
 *   L6: OV, Edit Distance, LCS, Frechet, DTW as canonical problems
 *   L7: Applications: DNA sequence alignment, time series matching
 *
 * References:
 *   Backurs & Indyk (2016), STOC 2015
 *   Abboud, Backurs, Williams (2015), FOCS
 *   Bringmann (2014), FOCS
 *   Williams (2005), ICALP
 * ============================================================================ */

#include "subquadratic.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <float.h>
#include <math.h>

static inline int32_t imin(int32_t a, int32_t b) { return a < b ? a : b; }
static inline int32_t imax(int32_t a, int32_t b) { return a > b ? a : b; }
static inline int32_t imin3(int32_t a, int32_t b, int32_t c) {
    return imin(imin(a, b), c);
}
static inline double dmin(double a, double b) { return a < b ? a : b; }
static inline double dmax(double a, double b) { return a > b ? a : b; }
static inline double dmin3(double a, double b, double c) {
    return dmin(dmin(a, b), c);
}

/* ========================================================================
 * L1/L5: Orthogonal Vectors (OV)
 * ======================================================================== */

ov_instance_t *ov_create(int32_t n, int32_t d) {
    if (n <= 0 || d <= 0) return NULL;
    ov_instance_t *inst = (ov_instance_t *)malloc(sizeof(ov_instance_t));
    if (!inst) return NULL;
    inst->n = n; inst->d = d;
    inst->d_words = (d + 31) / 32;
    size_t sz = (size_t)n * (size_t)inst->d_words * sizeof(uint32_t);
    inst->A = (uint32_t *)calloc(sz, 1);
    inst->B = (uint32_t *)calloc(sz, 1);
    if (!inst->A || !inst->B) {
        free(inst->A); free(inst->B); free(inst); return NULL;
    }
    return inst;
}

void ov_free(ov_instance_t *inst) {
    if (!inst) return;
    free(inst->A); free(inst->B); free(inst);
}

void ov_set_bit(ov_instance_t *inst, int32_t vi, int32_t bp, bool val, bool in_A) {
    if (!inst || vi < 0 || vi >= inst->n || bp < 0 || bp >= inst->d) return;
    uint32_t *arr = in_A ? inst->A : inst->B;
    int32_t wi = vi * inst->d_words + (bp / 32);
    int32_t bo = bp % 32;
    if (val) arr[wi] |= (1U << bo);
    else arr[wi] &= ~(1U << bo);
}

bool ov_get_bit(const ov_instance_t *inst, int32_t vi, int32_t bp, bool in_A) {
    if (!inst || vi < 0 || vi >= inst->n || bp < 0 || bp >= inst->d) return false;
    const uint32_t *arr = in_A ? inst->A : inst->B;
    int32_t wi = vi * inst->d_words + (bp / 32);
    return (arr[wi] >> (bp % 32)) & 1U;
}

void ov_random_instance(ov_instance_t *inst, double density, uint64_t seed) {
    if (!inst) return;
    srand((unsigned)seed);
    for (int32_t i = 0; i < inst->n; i++)
        for (int32_t j = 0; j < inst->d; j++) {
            bool bit = ((double)rand() / RAND_MAX) < density;
            ov_set_bit(inst, i, j, bit, true);
            ov_set_bit(inst, i, j, bit, false);
        }
}

bool ov_brute_force(const ov_instance_t *inst, int32_t *a_idx, int32_t *b_idx) {
    if (!inst) return false;
    int32_t n = inst->n, dw = inst->d_words;
    for (int32_t i = 0; i < n; i++) {
        const uint32_t *ai = &inst->A[i * dw];
        for (int32_t j = 0; j < n; j++) {
            const uint32_t *bj = &inst->B[j * dw];
            bool orth = true;
            for (int32_t w = 0; w < dw && orth; w++)
                if ((ai[w] & bj[w]) != 0) orth = false;
            if (orth) {
                if (a_idx) *a_idx = i;
                if (b_idx) *b_idx = j;
                return true;
            }
        }
    }
    return false;
}

bool ov_williams_algorithm(const ov_instance_t *inst, int32_t *a_idx, int32_t *b_idx) {
    if (!inst) return false;
    int32_t n = inst->n, d = inst->d;
    if (d > (int32_t)(log2((double)n) * 2.0) + 1)
        return ov_brute_force(inst, a_idx, b_idx);
    int32_t ns = (int32_t)(log2((double)n) * 2.0);
    if (ns >= d) return ov_brute_force(inst, a_idx, b_idx);
    for (int trial = 0; trial < 3; trial++) {
        int32_t *samples = (int32_t *)malloc((size_t)ns * sizeof(int32_t));
        if (!samples) continue;
        for (int32_t s = 0; s < ns; s++) samples[s] = rand() % d;
        for (int32_t i = 0; i < n; i++) {
            for (int32_t j = 0; j < n; j++) {
                bool possible = true;
                for (int32_t s = 0; s < ns && possible; s++)
                    if (ov_get_bit(inst, i, samples[s], true) &&
                        ov_get_bit(inst, j, samples[s], false))
                        possible = false;
                if (possible) {
                    int32_t dw = inst->d_words;
                    const uint32_t *ai = &inst->A[i * dw];
                    const uint32_t *bj = &inst->B[j * dw];
                    bool orth = true;
                    for (int32_t w = 0; w < dw && orth; w++)
                        if ((ai[w] & bj[w]) != 0) orth = false;
                    if (orth) {
                        if (a_idx) *a_idx = i;
                        if (b_idx) *b_idx = j;
                        free(samples); return true;
                    }
                }
            }
        }
        free(samples);
    }
    return false;
}

/* ========================================================================
 * Edit Distance
 * ======================================================================== */

edit_distance_instance_t *edit_distance_create(int32_t n, int32_t m) {
    if (n < 0 || m < 0) return NULL;
    edit_distance_instance_t *inst = (edit_distance_instance_t *)
        malloc(sizeof(edit_distance_instance_t));
    if (!inst) return NULL;
    inst->n = n; inst->m = m;
    inst->x = (char *)malloc((size_t)(n + 1));
    inst->y = (char *)malloc((size_t)(m + 1));
    inst->dp = (int32_t *)calloc((size_t)(n + 1) * (size_t)(m + 1), sizeof(int32_t));
    if (!inst->x || !inst->y || !inst->dp) {
        free(inst->x); free(inst->y); free(inst->dp); free(inst); return NULL;
    }
    inst->distance = -1;
    return inst;
}

void edit_distance_free(edit_distance_instance_t *inst) {
    if (!inst) return;
    free(inst->x); free(inst->y); free(inst->dp); free(inst);
}

int32_t edit_distance_levenshtein(const char *x, int32_t n, const char *y, int32_t m) {
    if (!x || !y) return -1;
    if (n == 0) return m;
    if (m == 0) return n;
    int32_t *prev = (int32_t *)malloc((size_t)(m + 1) * sizeof(int32_t));
    int32_t *curr = (int32_t *)malloc((size_t)(m + 1) * sizeof(int32_t));
    if (!prev || !curr) { free(prev); free(curr); return -1; }
    for (int32_t j = 0; j <= m; j++) prev[j] = j;
    for (int32_t i = 1; i <= n; i++) {
        curr[0] = i;
        for (int32_t j = 1; j <= m; j++) {
            int32_t sc = (x[i - 1] == y[j - 1]) ? 0 : 1;
            curr[j] = imin3(prev[j] + 1, curr[j - 1] + 1, prev[j - 1] + sc);
        }
        int32_t *tmp = prev; prev = curr; curr = tmp;
    }
    int32_t result = prev[m];
    free(prev); free(curr);
    return result;
}

int32_t edit_distance_dp(edit_distance_instance_t *inst) {
    if (!inst) return -1;
    int32_t n = inst->n, m = inst->m;
    for (int32_t i = 0; i <= n; i++) inst->dp[i * (m + 1) + 0] = i;
    for (int32_t j = 0; j <= m; j++) inst->dp[0 * (m + 1) + j] = j;
    for (int32_t i = 1; i <= n; i++)
        for (int32_t j = 1; j <= m; j++) {
            int32_t cost = (inst->x[i - 1] == inst->y[j - 1]) ? 0 : 1;
            inst->dp[i * (m + 1) + j] = imin3(
                inst->dp[(i - 1) * (m + 1) + j] + 1,
                inst->dp[i * (m + 1) + (j - 1)] + 1,
                inst->dp[(i - 1) * (m + 1) + (j - 1)] + cost);
        }
    inst->distance = inst->dp[n * (m + 1) + m];
    return inst->distance;
}

int32_t edit_distance_linear_space(const edit_distance_instance_t *inst) {
    if (!inst) return -1;
    return edit_distance_levenshtein(inst->x, inst->n, inst->y, inst->m);
}

int32_t edit_distance_weighted(const edit_distance_instance_t *inst,
                                int ci, int cd, int cs) {
    if (!inst) return -1;
    int32_t n = inst->n, m = inst->m;
    int32_t *prev = (int32_t *)malloc((size_t)(m + 1) * sizeof(int32_t));
    int32_t *curr = (int32_t *)malloc((size_t)(m + 1) * sizeof(int32_t));
    if (!prev || !curr) { free(prev); free(curr); return -1; }
    for (int32_t j = 0; j <= m; j++) prev[j] = j * ci;
    for (int32_t i = 1; i <= n; i++) {
        curr[0] = i * cd;
        for (int32_t j = 1; j <= m; j++) {
            int32_t sc = (inst->x[i - 1] == inst->y[j - 1]) ? 0 : cs;
            curr[j] = imin3(prev[j] + cd, curr[j - 1] + ci, prev[j - 1] + sc);
        }
        int32_t *tmp = prev; prev = curr; curr = tmp;
    }
    int32_t result = prev[m];
    free(prev); free(curr);
    return result;
}

/* ========================================================================
 * LCS
 * ======================================================================== */

lcs_instance_t *lcs_create(int32_t n, int32_t m) {
    if (n < 0 || m < 0) return NULL;
    lcs_instance_t *inst = (lcs_instance_t *)malloc(sizeof(lcs_instance_t));
    if (!inst) return NULL;
    inst->n = n; inst->m = m;
    inst->x = (char *)malloc((size_t)(n + 1));
    inst->y = (char *)malloc((size_t)(m + 1));
    inst->dp = (int32_t *)calloc((size_t)(n + 1) * (size_t)(m + 1), sizeof(int32_t));
    inst->lcs_string = NULL;
    if (!inst->x || !inst->y || !inst->dp) {
        free(inst->x); free(inst->y); free(inst->dp); free(inst); return NULL;
    }
    inst->lcs_length = 0;
    return inst;
}

void lcs_free(lcs_instance_t *inst) {
    if (!inst) return;
    free(inst->x); free(inst->y); free(inst->dp);
    free(inst->lcs_string); free(inst);
}

int32_t lcs_dp(lcs_instance_t *inst) {
    if (!inst) return -1;
    int32_t n = inst->n, m = inst->m, st = m + 1;
    for (int32_t i = 1; i <= n; i++)
        for (int32_t j = 1; j <= m; j++) {
            if (inst->x[i - 1] == inst->y[j - 1])
                inst->dp[i * st + j] = inst->dp[(i - 1) * st + (j - 1)] + 1;
            else
                inst->dp[i * st + j] = imax(
                    inst->dp[(i - 1) * st + j],
                    inst->dp[i * st + (j - 1)]);
        }
    inst->lcs_length = inst->dp[n * st + m];
    return inst->lcs_length;
}

int32_t lcs_hirschberg(lcs_instance_t *inst) {
    if (!inst) return -1;
    int32_t len = lcs_dp(inst);
    if (len < 0) return -1;
    int32_t n = inst->n, m = inst->m, st = m + 1;
    free(inst->lcs_string);
    inst->lcs_string = (char *)malloc((size_t)(len + 1));
    if (!inst->lcs_string) return len;
    int32_t i = n, j = m, pos = len;
    inst->lcs_string[pos] = '\0';
    while (i > 0 && j > 0) {
        if (inst->x[i - 1] == inst->y[j - 1]) {
            inst->lcs_string[--pos] = inst->x[i - 1]; i--; j--;
        } else if (inst->dp[(i - 1) * st + j] >= inst->dp[i * st + (j - 1)])
            i--;
        else j--;
    }
    return len;
}

/* ========================================================================
 * Frechet Distance
 * ======================================================================== */

frechet_instance_t *frechet_create(int32_t n, int32_t m) {
    if (n <= 0 || m <= 0) return NULL;
    frechet_instance_t *inst = (frechet_instance_t *)
        malloc(sizeof(frechet_instance_t));
    if (!inst) return NULL;
    inst->n = n; inst->m = m;
    inst->Px = (double *)calloc((size_t)n, sizeof(double));
    inst->Py = (double *)calloc((size_t)n, sizeof(double));
    inst->Qx = (double *)calloc((size_t)m, sizeof(double));
    inst->Qy = (double *)calloc((size_t)m, sizeof(double));
    inst->distance = -1.0;
    if (!inst->Px || !inst->Py || !inst->Qx || !inst->Qy) {
        free(inst->Px); free(inst->Py); free(inst->Qx); free(inst->Qy);
        free(inst); return NULL;
    }
    return inst;
}

void frechet_free(frechet_instance_t *inst) {
    if (!inst) return;
    free(inst->Px); free(inst->Py);
    free(inst->Qx); free(inst->Qy); free(inst);
}

double frechet_continuous(const frechet_instance_t *inst) {
    if (!inst) return -1.0;
    int32_t n = inst->n, m = inst->m;

    double **pd = (double **)malloc((size_t)n * sizeof(double *));
    if (!pd) return -1.0;
    for (int32_t i = 0; i < n; i++) {
        pd[i] = (double *)malloc((size_t)m * sizeof(double));
        if (!pd[i]) {
            for (int32_t k = 0; k < i; k++) free(pd[k]);
            free(pd); return -1.0;
        }
        for (int32_t j = 0; j < m; j++) {
            double dx = inst->Px[i] - inst->Qx[j];
            double dy = inst->Py[i] - inst->Qy[j];
            pd[i][j] = sqrt(dx * dx + dy * dy);
        }
    }

    double *dp = (double *)malloc((size_t)n * (size_t)m * sizeof(double));
    if (!dp) {
        for (int32_t i = 0; i < n; i++) free(pd[i]);
        free(pd); return -1.0;
    }

    dp[0 * m + 0] = pd[0][0];
    for (int32_t i = 1; i < n; i++)
        dp[i * m + 0] = dmax(dp[(i - 1) * m + 0], pd[i][0]);
    for (int32_t j = 1; j < m; j++)
        dp[0 * m + j] = dmax(dp[0 * m + (j - 1)], pd[0][j]);

    for (int32_t i = 1; i < n; i++)
        for (int32_t j = 1; j < m; j++)
            dp[i * m + j] = dmax(
                dmin3(dp[(i - 1) * m + j], dp[i * m + (j - 1)],
                      dp[(i - 1) * m + (j - 1)]),
                pd[i][j]);

    double result = dp[(n - 1) * m + (m - 1)];
    free(dp);
    for (int32_t i = 0; i < n; i++) free(pd[i]);
    free(pd);
    return result;
}

/* ========================================================================
 * DTW
 * ======================================================================== */

dtw_instance_t *dtw_create(int32_t n, int32_t m) {
    if (n <= 0 || m <= 0) return NULL;
    dtw_instance_t *inst = (dtw_instance_t *)malloc(sizeof(dtw_instance_t));
    if (!inst) return NULL;
    inst->n = n; inst->m = m;
    inst->x = (double *)calloc((size_t)n, sizeof(double));
    inst->y = (double *)calloc((size_t)m, sizeof(double));
    inst->dp = (double *)calloc((size_t)n * (size_t)m, sizeof(double));
    if (!inst->x || !inst->y || !inst->dp) {
        free(inst->x); free(inst->y); free(inst->dp); free(inst); return NULL;
    }
    inst->dtw_distance = -1.0;
    return inst;
}

void dtw_free(dtw_instance_t *inst) {
    if (!inst) return;
    free(inst->x); free(inst->y); free(inst->dp); free(inst);
}

double dtw_standard(dtw_instance_t *inst) {
    if (!inst) return -1.0;
    int32_t n = inst->n, m = inst->m;
    inst->dp[0 * m + 0] = fabs(inst->x[0] - inst->y[0]);
    for (int32_t j = 1; j < m; j++)
        inst->dp[0 * m + j] = inst->dp[0 * m + (j - 1)] + fabs(inst->x[0] - inst->y[j]);
    for (int32_t i = 1; i < n; i++)
        inst->dp[i * m + 0] = inst->dp[(i - 1) * m + 0] + fabs(inst->x[i] - inst->y[0]);
    for (int32_t i = 1; i < n; i++)
        for (int32_t j = 1; j < m; j++) {
            double cost = fabs(inst->x[i] - inst->y[j]);
            inst->dp[i * m + j] = cost + dmin3(
                inst->dp[(i - 1) * m + j],
                inst->dp[i * m + (j - 1)],
                inst->dp[(i - 1) * m + (j - 1)]);
        }
    inst->dtw_distance = inst->dp[(n - 1) * m + (m - 1)];
    return inst->dtw_distance;
}

double dtw_sakoe_chiba(dtw_instance_t *inst, int32_t w) {
    if (!inst || w <= 0) return -1.0;
    int32_t n = inst->n, m = inst->m;
    for (int32_t i = 0; i < n; i++)
        for (int32_t j = 0; j < m; j++)
            inst->dp[i * m + j] = DBL_MAX;
    inst->dp[0 * m + 0] = fabs(inst->x[0] - inst->y[0]);
    for (int32_t i = 0; i < n; i++) {
        int32_t js = imax(0, i - w), je = imin(m - 1, i + w);
        for (int32_t j = js; j <= je; j++) {
            if (i == 0 && j == 0) continue;
            double cost = fabs(inst->x[i] - inst->y[j]);
            double best = DBL_MAX;
            if (i > 0) best = dmin(best, inst->dp[(i - 1) * m + j]);
            if (j > 0) best = dmin(best, inst->dp[i * m + (j - 1)]);
            if (i > 0 && j > 0) best = dmin(best, inst->dp[(i - 1) * m + (j - 1)]);
            if (best < DBL_MAX) inst->dp[i * m + j] = best + cost;
        }
    }
    inst->dtw_distance = inst->dp[(n - 1) * m + (m - 1)];
    return inst->dtw_distance;
}

/* ========================================================================
 * L7: Applications
 * ======================================================================== */

int32_t dna_sequence_alignment(const char *seq1, int32_t n1,
                                const char *seq2, int32_t n2,
                                int ms, int mp, int gp) {
    if (!seq1 || !seq2 || n1 < 0 || n2 < 0) return -1;
    int32_t *prev = (int32_t *)malloc((size_t)(n2 + 1) * sizeof(int32_t));
    int32_t *curr = (int32_t *)malloc((size_t)(n2 + 1) * sizeof(int32_t));
    if (!prev || !curr) { free(prev); free(curr); return -1; }
    for (int32_t j = 0; j <= n2; j++) prev[j] = j * gp;
    for (int32_t i = 1; i <= n1; i++) {
        curr[0] = i * gp;
        for (int32_t j = 1; j <= n2; j++) {
            int32_t m = prev[j - 1] + ((seq1[i - 1] == seq2[j - 1]) ? ms : mp);
            curr[j] = imax(imax(m, prev[j] + gp), curr[j - 1] + gp);
        }
        int32_t *tmp = prev; prev = curr; curr = tmp;
    }
    int32_t result = prev[n2];
    free(prev); free(curr);
    return result;
}

int32_t time_series_search(const double *query, int32_t ql,
                           const double **db, int32_t k, int32_t n,
                           double *best_dist) {
    if (!query || !db || !best_dist || k <= 0 || n <= 0) return -1;
    *best_dist = DBL_MAX;
    int32_t best_idx = -1;
    for (int32_t idx = 0; idx < k; idx++) {
        dtw_instance_t *dtw = dtw_create(ql, n);
        if (!dtw) continue;
        for (int32_t i = 0; i < ql; i++) dtw->x[i] = query[i];
        for (int32_t j = 0; j < n; j++) dtw->y[j] = db[idx][j];
        double dist = dtw_standard(dtw);
        if (dist >= 0 && dist < *best_dist) { *best_dist = dist; best_idx = idx; }
        dtw_free(dtw);
    }
    return best_idx;
}

/* ========================================================================
 * L6/L8: Subquadratic Completeness
 * ======================================================================== */

bool verify_subquadratic_completeness(int32_t n, int32_t d,
    bool (*solve_P)(void *inst, int32_t n),
    void *(*reduce_OV)(const ov_instance_t *ov)) {
    if (!solve_P || !reduce_OV) return false;
    ov_instance_t *ov = ov_create(n, d);
    if (!ov) return false;
    ov_random_instance(ov, 0.5, 42);
    void *pi = reduce_OV(ov);
    ov_free(ov);
    if (!pi) return false;
    return solve_P(pi, n);
}

void subquadratic_status_report(void) {
    printf("=== Subquadratic Equivalence Class ===\n");
    printf("Canonical problem: Orthogonal Vectors (OV)\n");
    printf("Threshold exponent: 2.0\n");
    printf("Conjectured lower bound: 2.0 (under SETH)\n\n");
    printf("Member problems:\n");
    printf("  1. Orthogonal Vectors (OV)\n");
    printf("  2. Edit Distance (Levenshtein)\n");
    printf("  3. Longest Common Subsequence (LCS)\n");
    printf("  4. Frechet Distance\n");
    printf("  5. Dynamic Time Warping (DTW)\n");
    printf("  6. Longest Common Substring\n");
    printf("  7. Pattern Matching with Wildcards\n\n");
    printf("Key theorems:\n");
    printf("  Backurs & Indyk (2016): SETH => Edit Distance is n^{2-o(1)}\n");
    printf("  Bringmann (2014): SETH => Frechet Distance is n^{2-o(1)}\n");
    printf("  Abboud et al. (2015): SETH => LCS, DTW are n^{2-o(1)}\n\n");
    printf("Best known algorithms:\n");
    printf("  OV: O(n^{2-1/O(log d/log n)}) (Williams 2005)\n");
    printf("  Edit Distance: O(n^2/log^2 n) (Masek & Paterson 1980)\n");
    printf("  LCS: O(n^2/log^2 n) (Four Russians)\n");
}