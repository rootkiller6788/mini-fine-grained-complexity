#ifndef FRECHET_DISTANCE_H
#define FRECHET_DISTANCE_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <math.h>

/* ============================================================================
 * mini-edit-distance-lcs-hardness: Frechet Distance
 *
 * Frechet distance (dog-leash distance) between two curves.
 * Conditional lower bound: O(n^{2-ε}) under SETH.
 *
 * Reference: Bringmann & Kunnemann (FOCS 2015)
 *            "Quadratic Conditional Lower Bounds for String Problems
 *            and Dynamic Time Warping"
 *
 * L1: Frechet distance definition
 * L2: Quadratic-time lower bound under SETH
 * L5: Discrete Frechet distance computation
 * ============================================================================ */

/* A point in 2D space */
typedef struct {
    double x;
    double y;
} point2d_t;

/* A polygonal curve: sequence of 2D points */
typedef struct {
    point2d_t *points;
    int32_t    num_points;
} curve_t;

/* ---- L1: Curve Operations ---- */

curve_t *curve_create(int32_t num_points);
void curve_destroy(curve_t *curve);
curve_t *curve_from_arrays(const double *x, const double *y, int32_t n);
double point_distance(const point2d_t *a, const point2d_t *b);
double curve_length(const curve_t *curve);

/* ---- L2: Frechet Distance ---- */

/* Discrete Frechet distance between two curves.
 * O(n*m) time, O(n*m) space using DP.
 * This is the dog-leash distance: the minimum leash length
 * needed for a person and dog to traverse their respective
 * curves while staying connected. */
double discrete_frechet_distance(const curve_t *p, const curve_t *q);

/* Continuous Frechet distance approximation via discretization.
 * Discretizes each curve segment into k sub-segments.
 * As k → ∞, approaches the true continuous Frechet distance. */
double continuous_frechet_approximation(const curve_t *p, const curve_t *q,
                                         int32_t discretization);

/* Frechet distance decision problem: is d_F(P,Q) ≤ δ?
 * Uses free-space diagram reachability.
 * O(n*m) time. */
bool frechet_decision(const curve_t *p, const curve_t *q, double delta);

/* ---- Dynamic Time Warping (DTW) ---- */

/* Dynamic Time Warping distance between two curves.
 * Similar to Frechet but allows one-to-many matching.
 * O(n*m) time, O(n*m) space. */
double dtw_distance(const curve_t *p, const curve_t *q);

/* DTW with Sakoe-Chiba band constraint (window size w).
 * O(w * min(n,m)) time. */
double dtw_banded(const curve_t *p, const curve_t *q, int32_t w);

/* ---- L4: Conditional Lower Bound ---- */

/* Verify Bringmann-Kunnemann reduction:
 * SAT(n) → OV(N,d) → Frechet(N')
 * If Frechet can be computed in O(N^{2-ε}) time,
 * then SAT can be solved in sub-exponential time. */
bool frechet_seth_hardness_verify(int32_t n, double epsilon);

/* Print Frechet distance SETH status report */
void print_frechet_seth_status(void);

#endif /* FRECHET_DISTANCE_H */
