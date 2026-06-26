#ifndef EQUIVALENCE_H
#define EQUIVALENCE_H
#include <stdint.h>
#include <stdbool.h>
#ifdef __cplusplus
extern "C" {
#endif
typedef enum { EQ_APSP=0,EQ_MIN_PLUS_PRODUCT=1,EQ_NEGATIVE_TRIANGLE=2,EQ_RADIUS=3,EQ_MEDIAN=4,EQ_BETWEENNESS=5,EQ_DIAMETER=6,EQ_3SUM=7,EQ_COLLINEAR=8,EQ_MIN_AREA_TRIANGLE=9,EQ_COUNT=10 } eq_problem_t;
typedef struct { eq_problem_t problem; const char *name; double trivial_bound; double conjectured_bound; double best_known; int32_t equivalence_class_id; } eq_problem_info_t;
const char *eq_problem_name(eq_problem_t p);
eq_problem_info_t eq_get_info(eq_problem_t p);
bool eq_are_equivalent(eq_problem_t a, eq_problem_t b);
void eq_print_equivalence_map(void);
bool eq_verify_hexagon(void);
#ifdef __cplusplus
}
#endif
#endif

/*
 * The subcubic equivalence hexagon (Williams & Williams 2010):
 *
 *        APSP ----------> Min-Plus Product
 *         |                     |
 *         v                     v
 *       Radius ----------> Negative Triangle
 *         |                     |
 *         v                     v
 *       Median --------> Betweenness Centrality
 *
 * Each directed edge represents a subcubic reduction.
 * The graph is strongly connected (all problems equivalent).
 */
/*
 * Fine-Grained Complexity Classes:
 *
 * Class 0 (Subcubic): APSP, Min-Plus Product, Negative Triangle,
 *                     Radius, Median, Betweenness, Diameter
 *   Barrier: O(n^{3-epsilon})
 *
 * Class 1 (Subquadratic): 3SUM, Collinearity, Min-Area Triangle
 *   Barrier: O(n^{2-epsilon})
 *
 * Bridge (VW 2009): Class 0 and Class 1 are subcubic-equivalent.
 */