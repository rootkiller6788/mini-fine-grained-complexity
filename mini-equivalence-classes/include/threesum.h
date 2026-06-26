#ifndef THREESUM_H
#define THREESUM_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <math.h>

/*
 * threesum.h -- 3SUM Equivalence Class
 *
 * The 3SUM equivalence class contains problems whose best known
 * algorithms run in O(n^2) time and for which no O(n^{2-epsilon})
 * algorithm is known for any epsilon > 0, in the decision tree model
 * or in general.
 *
 * Canonical problem: 3SUM
 * Given a set S of n integers, are there a,b,c in S with a + b + c = 0?
 *
 * 3SUM Conjecture: No O(n^{2-epsilon}) algorithm for 3SUM.
 *
 * The 3SUM conjecture is believed even in the integer regime with
 * values in [-n^3, n^3], but weaker than SETH-based conjectures in
 * the sense that 3SUM can be solved in O(n^2 / log^2 n) (Baran,
 * Demaine, Patrascu 2005) while subquadratic-equivalent problems
 * like OV have stronger SETH-based lower bounds.
 *
 * References:
 *   Gajentaan & Overmars (1995), "On a Class of O(n^2) Problems in
 *     Computational Geometry", Computational Geometry
 *   Baran, Demaine, Patrascu (2005), "Subquadratic Algorithms for 3SUM",
 *     WADS
 *   Gronlund & Pettie (2014), "Threesomes, Degenerates, and Love Triangles",
 *     FOCS (solves 3SUM in O(n^2 / (log n / log log n)^2/3) deterministic)
 *
 * L1: 3SUM problem definition, variants (3SUM', 3SUM+)
 * L2: 3SUM-hardness concept, 3SUM-equivalence
 * L3: Integer sets, geometric point configurations
 * L4: 3SUM conjecture, reductions to geometric problems
 * L5: Sorting-based 3SUM (O(n^2)), Gronlund-Pettie (O(n^2 / polylog n))
 * L6: 3SUM, 3SUM', Collinearity, Polygon Fitting
 */

/* ---- L1: 3SUM Instance ---- */

typedef struct {
    int32_t     n;             /* Number of integers */
    int64_t    *values;        /* Array of n integers */
    bool        found;         /* True if 3SUM solution found */
    int32_t     idx_a;         /* Index of first element (a) */
    int32_t     idx_b;         /* Index of second element (b) */
    int32_t     idx_c;         /* Index of third element (c) */
} threesum_instance_t;

/* ---- L1: 3SUM' (3SUM-prime) Instance ---- */

typedef struct {
    int32_t     n;             /* Number of integers in each of 3 sets */
    int64_t    *A;             /* First set, length n */
    int64_t    *B;             /* Second set, length n */
    int64_t    *C;             /* Third set, length n */
    bool        found;
    int32_t     ia;            /* Index in A */
    int32_t     ib;            /* Index in B */
    int32_t     ic;            /* Index in C */
} threesum_prime_instance_t;

/* ---- L1: Collinearity Instance (3SUM-hard) ---- */

typedef struct {
    int32_t     n;             /* Number of points */
    double     *x;             /* X coordinates, length n */
    double     *y;             /* Y coordinates, length n */
    bool        found;         /* True if 3 collinear points found */
    int32_t     i1, i2, i3;    /* Indices of collinear points */
} collinearity_instance_t;

/* ---- L2/L5: 3SUM Algorithms ---- */

threesum_instance_t *threesum_create(int32_t n);
void threesum_free(threesum_instance_t *inst);
void threesum_set(threesum_instance_t *inst, int32_t i, int64_t value);
int64_t threesum_get(const threesum_instance_t *inst, int32_t i);

/*
 * 3SUM via sorting + two-pointer technique.
 * First sort in O(n log n), then for each a, find b,c with b+c = -a
 * using two pointers in O(n). Total: O(n^2).
 *
 * This is the classic O(n^2) algorithm and the baseline for the
 * 3SUM conjecture.
 */
bool threesum_quadratic(threesum_instance_t *inst);

/*
 * 3SUM' (three distinct sets A, B, C).
 * Find a in A, b in B, c in C such that a + b + c = 0.
 * Reduces to standard 3SUM by S = A U (B shifted) U (C shifted).
 */
bool threesum_prime_quadratic(threesum_prime_instance_t *inst);

/*
 * Gronlund-Pettie (2014) deterministic 3SUM.
 * Time: O(n^2 * (log log n / log n)^(2/3)).
 * Currently the best known deterministic algorithm for 3SUM.
 *
 * This is a simplified version that demonstrates the key ideas.
 */
bool threesum_gronlund_pettie(threesum_instance_t *inst);

/*
 * 3SUM with small integer range [0, M].
 * When M = O(n^2), can be solved in O(n log n) using FFT.
 * When M = O(n^3), the 3SUM conjecture is believed to hold.
 */
bool threesum_bounded(threesum_instance_t *inst, int64_t M);

/* ---- L6: 3SUM-hard problems ---- */

collinearity_instance_t *collinearity_create(int32_t n);
void collinearity_free(collinearity_instance_t *inst);

/*
 * Collinearity detection: do 3 points lie on a common line?
 * This is 3SUM-hard (Gajentaan & Overmars 1995).
 * Uses the standard reduction: points (x,y) -> x values, with
 * slope a -> integer encoding.
 *
 * Time: O(n^2) using 3SUM.
 */
bool collinearity_via_3sum(collinearity_instance_t *inst);

/*
 * Set disjointness of three sets A, B, C of real numbers.
 * Determine if there exist a in A, b in B, c in C all equal.
 * This is equivalent to 3SUM: a + b = c maps to a' + b' + c' = 0.
 */
bool three_set_equality(const int64_t *A, const int64_t *B, const int64_t *C, int32_t n);

/* ---- L7: Applications ---- */

/*
 * Computational geometry: point location with 3SUM-hard preprocessing.
 * Detecting if three points form an axis-aligned right triangle.
 */
bool right_triangle_detection(const double *x, const double *y, int32_t n,
                               int32_t *i1, int32_t *i2, int32_t *i3);

/*
 * Motion planning: detecting if three obstacles are collinear
 * (important for visibility graph construction).
 *
 * NASA space robotics: path planning for Mars rover obstacle avoidance.
 * Uses collinearity detection in O(n^2) via 3SUM reduction.
 */
bool mars_rover_collinearity_obstacle_check(
    const double *obs_x, const double *obs_y, int32_t num_obstacles);

/* ---- L8: Advanced ---- */

/*
 * Verify that a problem P is 3SUM-complete:
 * 1. P can be solved in O(n^2) time.
 * 2. 3SUM reduces to P in O(n^{2-epsilon}) time.
 */
bool verify_3sum_completeness(int32_t n,
    bool (*solve_P)(void *instance, int32_t n),
    void *(*reduce_3SUM_to_P)(const threesum_instance_t *ts));

/*
 * Print a summary of the 3SUM equivalence class.
 */
void threesum_status_report(void);

#endif /* THREESUM_H */