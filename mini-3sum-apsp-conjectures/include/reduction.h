#ifndef REDUCTION_H
#define REDUCTION_H
#include <stdint.h>
#include <stdbool.h>
#include "threesum.h"
#include "apsp.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef enum { RED_3SUM_TO_APSP, RED_APSP_TO_3SUM, RED_NEG_TRIANGLE, RED_MIN_PLUS, RED_COUNT } reduction_type_t;
typedef struct { reduction_type_t type; const char *name; const char *source_problem; const char *target_problem; double time_overhead; const char *description; } reduction_info_t;
typedef struct { reduction_type_t type; int32_t n_original; int32_t n_reduced; double *reduced_data; void *auxiliary; } reduction_result_t;
reduction_info_t reduction_get_info(reduction_type_t t);
reduction_result_t *reduction_3sum_to_apsp(const ts_instance_t *ts);
reduction_result_t *reduction_apsp_to_3sum(const apsp_graph_t *g);
bool reduction_verify_subcubic_equivalence(void);
void reduction_result_destroy(reduction_result_t *r);
void reduction_print_info(reduction_type_t t);
void reduction_print_chain(void);
void reduction_print_complexity_landscape(void);
void reduction_print_fine_grained_map(void);
void reduction_self_test(void);
void reduction_verify_all_equivalences(void);
void reduction_print_vw2009_theorem(void);
bool reduction_verify_3sum_to_apsp_direct(const ts_instance_t *inst);
bool reduction_is_subcubic_preserving(reduction_type_t t);
double reduction_asymptotic_overhead(reduction_type_t t, int32_t n);
double reduction_time_overhead_analysis(int32_t n, reduction_type_t t);
ts_instance_t *reduction_apsp_hard_to_3sum(apsp_hardness_class_t cls, int32_t n);
bool reduction_3sum_hard_to_apsp(ts_hardness_class_t cls, const ts_instance_t *inst, apsp_graph_t **out);
void reduction_explain_fine_grained(void);
#ifdef __cplusplus
}
#endif
#endif

/*
 * ============================================================
 * Reductions: Extended Documentation
 * ============================================================
 *
 * Fine-grained reductions differ from classical polynomial-time
 * reductions in several key ways:
 *
 * 1. Precision: A fine-grained reduction must preserve not just
 *    polynomial-time solvability, but the EXACT exponent in the
 *    running time (up to o(1) in the exponent).
 *
 * 2. Overhead: The reduction itself must run in time strictly
 *    less than the conjectured lower bound. For 3SUM/APSP:
 *    - 3SUM -> APSP: O(n) reduction (well below O(n^2))
 *    - APSP -> 3SUM: O(n^2) reduction (well below O(n^3))
 *
 * 3. Gap amplification: A reduction that blows up the input
 *    by too much may not preserve subcubic equivalence.
 *
 * Subcubic reduction definition (VW 2009):
 *   A reduction from problem A(n) to problem B(m) is
 *   "subcubic" if:
 *   (a) m = O(n) (linear blowup), and
 *   (b) The reduction runs in O(n^{3-delta}) time for some
 *       delta > 0.
 *
 *   Then: B in O(m^{3-epsilon}) implies A in O(n^{3-epsilon}).
 *
 * Key reductions in the subcubic equivalence web:
 *   3SUM <-> APSP <-> Min-Plus Product <-> Negative Triangle
 *   <-> Radius <-> Median <-> Betweenness Centrality
 *
 * Each arrow represents an O(n^{3-delta})-time reduction for
 * some delta > 0, establishing that truly subcubic algorithms
 * for any problem imply truly subcubic algorithms for all.
 *
 * For further reading:
 *   - Vassilevska-Williams & Williams (2010): FOCS 2010
 *   - Bringmann (2019): "Fine-Grained Complexity Tutorial"
 */
/*
 * The Fine-Grained Reduction Web:
 *
 *               3SUM --(O(n))--> APSP
 *                 ^                |
 *                 |                v
 *            Collinearity    Min-Plus Product
 *                 ^                |
 *                 |                v
 *          Min Area Tri    Negative Triangle
 *                                  |
 *                                  v
 *                             Graph Radius
 *                                  |
 *                                  v
 *                             Graph Median
 *
 * All arrows are subcubic reductions. The web is strongly
 * connected: solving ANY problem faster would solve ALL faster.
 */