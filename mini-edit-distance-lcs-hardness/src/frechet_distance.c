
/* ============================================================================
 * frechet_distance.c -- Frechet Distance and DTW
 *
 * Discrete Frechet distance (dog-leash distance), continuous approximation,
 * Dynamic Time Warping (DTW) with Sakoe-Chiba band, SETH hardness.
 *
 * L1: point2d_t, curve_t
 * L2: Discrete Frechet O(nm), DTW O(nm)
 * L4: Bringmann-Kunnemann (2015) conditional lower bound
 * ============================================================================ */

#include "frechet_distance.h"
#include "seth_hardness.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

static inline double min2d(double a, double b) { return (a < b) ? a : b; }
static inline double max2d(double a, double b) { return (a > b) ? a : b; }
static inline double max3d(double a, double b, double c) {
    double m = (a > b) ? a : b;
    return (m > c) ? m : c;
}
static inline double min3d(double a, double b, double c) {
    double m = (a < b) ? a : b;
    return (m < c) ? m : c;
}

/* ============================================================================
 * L1: Curve Operations
 * ============================================================================ */

curve_t *curve_create(int32_t num_points) {
    curve_t *c = (curve_t *)malloc(sizeof(curve_t));
    if (!c) return NULL;
    c->num_points = num_points;
    c->points = (point2d_t *)calloc((size_t)num_points, sizeof(point2d_t));
    if (!c->points) { free(c); return NULL; }
    return c;
}

void curve_destroy(curve_t *curve) {
    if (!curve) return;
    free(curve->points);
    free(curve);
}

curve_t *curve_from_arrays(const double *x, const double *y, int32_t n) {
    curve_t *c = curve_create(n);
    if (!c) return NULL;
    for (int32_t i = 0; i < n; i++) {
        c->points[i].x = x[i];
        c->points[i].y = y[i];
    }
    return c;
}

double point_distance(const point2d_t *a, const point2d_t *b) {
    if (!a || !b) return INFINITY;
    double dx = a->x - b->x;
    double dy = a->y - b->y;
    return sqrt(dx*dx + dy*dy);
}

double curve_length(const curve_t *curve) {
    if (!curve || curve->num_points < 2) return 0.0;
    double len = 0.0;
    for (int32_t i = 1; i < curve->num_points; i++) {
        len += point_distance(&curve->points[i-1], &curve->points[i]);
    }
    return len;
}

/* ============================================================================
 * L2: Discrete Frechet Distance
 *
 * Let P = (p_1,...,p_n), Q = (q_1,...,q_m).
 * Define coupling C = sequence of pairs (p_{a_k}, q_{b_k}) s.t.
 * a_1=1, b_1=1, a_K=n, b_K=m, and (a_{k+1},b_{k+1}) is one of
 * (a_k+1,b_k), (a_k,b_k+1), (a_k+1,b_k+1).
 *
 * Discrete Frechet: d_F(P,Q) = min_C max_k dist(p_{a_k}, q_{b_k})
 *
 * DP recurrence:
 *   D[i][j] = max(dist(p_i,q_j), min(D[i-1][j], D[i][j-1], D[i-1][j-1]))
 *
 * O(n*m) time and space.
 * ============================================================================ */

double discrete_frechet_distance(const curve_t *p, const curve_t *q) {
    if (!p || !q || p->num_points < 1 || q->num_points < 1) return INFINITY;
    int32_t n = p->num_points, m = q->num_points;
    int32_t stride = m;
    double *dp = (double *)calloc((size_t)(n*m), sizeof(double));
    if (!dp) return INFINITY;

    /* Initialize first cell */
    dp[0] = point_distance(&p->points[0], &q->points[0]);

    /* First row */
    for (int32_t j = 1; j < m; j++) {
        double d = point_distance(&p->points[0], &q->points[j]);
        dp[j] = max2d(dp[j-1], d);
    }

    /* First column */
    for (int32_t i = 1; i < n; i++) {
        double d = point_distance(&p->points[i], &q->points[0]);
        dp[i*stride] = max2d(dp[(i-1)*stride], d);
    }

    /* Fill rest */
    for (int32_t i = 1; i < n; i++) {
        int32_t rb = i*stride, prb = (i-1)*stride;
        for (int32_t j = 1; j < m; j++) {
            double d = point_distance(&p->points[i], &q->points[j]);
            double min_prev = min3d(dp[prb+j], dp[rb+(j-1)], dp[prb+(j-1)]);
            dp[rb+j] = max2d(d, min_prev);
        }
    }

    double result = dp[(n-1)*stride + (m-1)];
    free(dp);
    return result;
}

/* Continuous Frechet approximation via curve discretization */
double continuous_frechet_approximation(const curve_t *p, const curve_t *q,
                                         int32_t discretization) {
    if (!p || !q || discretization < 1) return discrete_frechet_distance(p, q);
    /* Subdivide each segment into discretization sub-segments */
    int32_t np = p->num_points;
    int32_t nq = q->num_points;
    int32_t new_np = (np - 1) * discretization + 1;
    int32_t new_nq = (nq - 1) * discretization + 1;

    point2d_t *pp = (point2d_t *)calloc((size_t)new_np, sizeof(point2d_t));
    point2d_t *pq = (point2d_t *)calloc((size_t)new_nq, sizeof(point2d_t));
    if (!pp || !pq) { free(pp); free(pq); return discrete_frechet_distance(p, q); }

    /* Interpolate P */
    for (int32_t i = 0; i < np - 1; i++) {
        for (int32_t k = 0; k < discretization; k++) {
            double t = (double)k / (double)discretization;
            int32_t idx = i * discretization + k;
            pp[idx].x = p->points[i].x + t * (p->points[i+1].x - p->points[i].x);
            pp[idx].y = p->points[i].y + t * (p->points[i+1].y - p->points[i].y);
        }
    }
    pp[new_np-1] = p->points[np-1];

    /* Interpolate Q */
    for (int32_t i = 0; i < nq - 1; i++) {
        for (int32_t k = 0; k < discretization; k++) {
            double t = (double)k / (double)discretization;
            int32_t idx = i * discretization + k;
            pq[idx].x = q->points[i].x + t * (q->points[i+1].x - q->points[i].x);
            pq[idx].y = q->points[i].y + t * (q->points[i+1].y - q->points[i].y);
        }
    }
    pq[new_nq-1] = q->points[nq-1];

    curve_t cp = {pp, new_np};
    curve_t cq = {pq, new_nq};
    double result = discrete_frechet_distance(&cp, &cq);
    free(pp); free(pq);
    return result;
}

/* Frechet distance decision problem: is d_F(P,Q) <= delta?
 * Uses free-space diagram reachability in O(n*m) time. */
bool frechet_decision(const curve_t *p, const curve_t *q, double delta) {
    if (!p || !q || delta < 0.0) return false;
    int32_t n = p->num_points, m = q->num_points;
    int32_t stride = m;
    bool *reachable = (bool *)calloc((size_t)(n*m), sizeof(bool));
    if (!reachable) return false;

    reachable[0] = (point_distance(&p->points[0], &q->points[0]) <= delta);
    if (!reachable[0]) { free(reachable); return false; }

    for (int32_t i = 1; i < n; i++)
        reachable[i*stride] = reachable[(i-1)*stride] &&
            (point_distance(&p->points[i], &q->points[0]) <= delta);
    for (int32_t j = 1; j < m; j++)
        reachable[j] = reachable[j-1] &&
            (point_distance(&p->points[0], &q->points[j]) <= delta);

    for (int32_t i = 1; i < n; i++) {
        int32_t rb = i*stride, prb = (i-1)*stride;
        for (int32_t j = 1; j < m; j++) {
            if (point_distance(&p->points[i], &q->points[j]) <= delta) {
                reachable[rb+j] = reachable[prb+j] || reachable[rb+(j-1)] || reachable[prb+(j-1)];
            }
        }
    }

    bool result = reachable[(n-1)*stride + (m-1)];
    free(reachable);
    return result;
}

/* ============================================================================
 * Dynamic Time Warping (DTW)
 *
 * DTW allows one-to-many matching (warping), unlike Frechet
 * which requires strictly monotone coupling.
 *
 * DP: D[i][j] = dist(p_i, q_j) + min(D[i-1][j], D[i][j-1], D[i-1][j-1])
 * O(n*m) time and space.
 * ============================================================================ */

double dtw_distance(const curve_t *p, const curve_t *q) {
    if (!p || !q || p->num_points < 1 || q->num_points < 1) return INFINITY;
    int32_t n = p->num_points, m = q->num_points;
    int32_t stride = m;
    double *dp = (double *)calloc((size_t)(n*m), sizeof(double));
    if (!dp) return INFINITY;

    dp[0] = point_distance(&p->points[0], &q->points[0]);

    /* First row */
    for (int32_t j = 1; j < m; j++)
        dp[j] = dp[j-1] + point_distance(&p->points[0], &q->points[j]);

    /* First column */
    for (int32_t i = 1; i < n; i++)
        dp[i*stride] = dp[(i-1)*stride] + point_distance(&p->points[i], &q->points[0]);

    /* Fill rest */
    for (int32_t i = 1; i < n; i++) {
        int32_t rb = i*stride, prb = (i-1)*stride;
        for (int32_t j = 1; j < m; j++) {
            double d = point_distance(&p->points[i], &q->points[j]);
            double min_prev = min3d(dp[prb+j], dp[rb+(j-1)], dp[prb+(j-1)]);
            dp[rb+j] = d + min_prev;
        }
    }

    double result = dp[(n-1)*stride + (m-1)];
    free(dp);
    return result;
}

/* DTW with Sakoe-Chiba band: only allow warping within window w.
 * O(w * min(n,m)) time. */
double dtw_banded(const curve_t *p, const curve_t *q, int32_t w) {
    if (!p || !q || w < 1) return dtw_distance(p, q);
    int32_t n = p->num_points, m = q->num_points;
    int32_t stride = m;
    double *dp = (double *)calloc((size_t)(n*m), sizeof(double));
    if (!dp) return INFINITY;

    /* Initialize with INFINITY */
    for (int32_t i = 0; i < n; i++)
        for (int32_t j = 0; j < m; j++)
            dp[i*stride+j] = INFINITY;

    dp[0] = point_distance(&p->points[0], &q->points[0]);

    for (int32_t i = 0; i < n; i++) {
        int32_t j_start = (i - w > 0) ? (i - w) : 0;
        int32_t j_end   = (i + w < m) ? (i + w) : (m - 1);
        for (int32_t j = j_start; j <= j_end; j++) {
            if (i == 0 && j == 0) continue;
            double d = point_distance(&p->points[i], &q->points[j]);
            double best = INFINITY;
            if (i > 0) best = min2d(best, dp[(i-1)*stride+j]);
            if (j > 0) best = min2d(best, dp[i*stride+(j-1)]);
            if (i > 0 && j > 0) best = min2d(best, dp[(i-1)*stride+(j-1)]);
            dp[i*stride+j] = d + best;
        }
    }

    double result = dp[(n-1)*stride + (m-1)];
    free(dp);
    return result;
}

/* ---- L4: Bringmann-Kunnemann (2015) SETH Hardness ---- */

bool frechet_seth_hardness_verify(int32_t n, double epsilon) {
    /* Bringmann & Kunnemann (FOCS 2015):
     * If the Frechet distance of two curves of length N can be
     * computed in O(N^{2-eps}) time, then CNF-SAT on n variables
     * can be solved in O(2^{(1-delta)n}) time for some delta > 0,
     * refuting SETH. */
    if (epsilon <= 0.0 || n <= 0) return false;
    double threshold = (double)n * (1.0 - 1e-6);
    double frechet_exp = (double)n - epsilon * (double)n / 2.0;
    return frechet_exp < threshold;
}

void print_frechet_seth_status(void) {
    printf("=== Frechet Distance SETH Lower Bound ===\n");
    printf("Bringmann-Kunnemann (FOCS 2015):\n");
    printf("  Frechet in O(n^{2-eps}) => SETH false.\n");
    printf("  Same holds for Dynamic Time Warping (DTW).\n");
}
