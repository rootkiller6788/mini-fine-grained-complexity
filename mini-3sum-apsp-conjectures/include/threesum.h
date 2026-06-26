#ifndef THREESUM_H
#define THREESUM_H
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#ifdef __cplusplus
extern "C" {
#endif
typedef int64_t ts_elem_t;
typedef struct { ts_elem_t *vals; int32_t n; } ts_instance_t;
typedef struct { int32_t i, j, k; } ts_triple_t;
typedef struct { bool found; int32_t count; ts_triple_t *triples; int32_t capacity; } ts_result_t;
typedef enum { TS3_HARD_BASELINE=0, TS3_HARD_COLLINEAR=1, TS3_HARD_POLYGON=2, TS3_HARD_SEPARATOR=3, TS3_HARD_INTERSECTION=4, TS3_HARD_MOTION=5, TS3_HARD_COUNT=6 } ts_hardness_class_t;
typedef struct { const char *statement; double epsilon_min; int32_t model_type; const char *history; } ts_conjecture_t;
typedef struct { int64_t pairs_checked; int64_t comparisons; int64_t hash_probes; } ts_triple_space_t;
typedef struct { ts_elem_t *keys; int32_t *positions; bool *occupied; int32_t capacity; int32_t size; int64_t probes; } ts_hash_table_t;
typedef struct { double x; double y; } ts_point_t;
typedef struct { int64_t dx; int64_t dy; } ts_slope_t;
ts_result_t ts_naive_O_n3(const ts_instance_t *inst);
ts_result_t ts_quadratic_hash(const ts_instance_t *inst);
ts_result_t ts_quadratic_sort_bsearch(ts_instance_t *inst);
ts_result_t ts_quadratic_two_pointer(ts_instance_t *inst);
ts_result_t ts_discrete_3sum(const ts_elem_t *A, int32_t na, const ts_elem_t *B, int32_t nb, const ts_elem_t *C, int32_t nc);
bool ts_collinearity_test(const double *xs, const double *ys, int32_t n);
ts_conjecture_t ts_conjecture_status(void);
ts_result_t ts_gronlund_pettie_2014(ts_instance_t *inst);
typedef struct { const char *problem_name; double (*input_size_map)(int32_t n); double (*time_bound)(int32_t n); bool (*correctness)(void); } ts_fine_grained_reduction_t;
ts_instance_t *ts_instance_create(const ts_elem_t *vals, int32_t n);
ts_instance_t *ts_instance_clone(const ts_instance_t *inst);
void ts_instance_destroy(ts_instance_t *inst);
void ts_result_init(ts_result_t *res);
void ts_result_free(ts_result_t *res);
void ts_result_add_triple(ts_result_t *res, int32_t i, int32_t j, int32_t k);
void ts_print_instance(const ts_instance_t *inst);
void ts_print_result(const ts_result_t *res);
const char *ts_hardness_class_name(ts_hardness_class_t cls);
double ts_min_area_triangle(const double *xs, const double *ys, int32_t n, int32_t *out_i, int32_t *out_j, int32_t *out_k);
int64_t ts_count_zero_sum_triples(ts_instance_t *inst);
bool ts_3sum_hard_problem_verify(ts_hardness_class_t cls, const double *xs, const double *ys, int32_t n);
bool ts_is_3sum_hard(int32_t n, double observed_exponent);
ts_instance_t *ts_generate_random_instance(int32_t n, int64_t range, uint64_t seed);
#ifdef __cplusplus
}
#endif
#endif

/*
 * ============================================================
 * Extended Documentation
 * ============================================================
 *
 * The 3SUM problem has a rich history in computational geometry
 * and fine-grained complexity:
 *
 * 1977: Gajentaan and Overmars survey O(n^2) problems in
 *       computational geometry, noticing that many problems
 *       (collinearity, min-area triangle, etc.) seem to require
 *       quadratic time. They conjecture that 3SUM itself requires
 *       Omega(n^2) time in the algebraic decision tree model.
 *
 * 1995: Erickson proves an Omega(n^2) lower bound for 3SUM in
 *       the 3-linear decision tree model (each query is
 *       a*x_i + b*x_j + c*x_k ? 0). This is a restricted model,
 *       but it provides evidence for the conjecture.
 *
 * 2005: Ailon and Chazelle extend the lower bound to the
 *       k-linear decision tree model for any constant k.
 *
 * 2014: Gronlund and Pettie (FOCS 2014) give the first algorithm
 *       that beats n^2 by a super-polylogarithmic factor:
 *       O(n^2 * (log log n / log n)^{2/3}).
 *       Key technique: Fredman-Komlos-Szemeredi hashing plus
 *       dominance counting in 2D.
 *
 * 2017: Gold and Sharir (STOC 2017) improve to
 *       O(n^2 * log log n / log n).
 *       Key technique: improved polynomial partitioning.
 *
 * Despite these improvements, no algorithm achieves O(n^{1.999})
 * for any fixed epsilon > 0. The 3SUM Conjecture, in its modern
 * form, states that 3SUM requires n^{2-o(1)} time.
 *
 * The 3SUM Conjecture has implications beyond computational
 * geometry. Under subcubic reductions (Vassilevska-Williams 2009),
 * it is equivalent to the APSP Conjecture: 3SUM is in
 * O(n^{2-epsilon}) iff APSP is in O(n^{3-epsilon}).
 *
 * The reduction 3SUM -> APSP works as follows:
 *   Given numbers x_1, ..., x_n, construct a tripartite graph
 *   with 3n vertices (A_1..A_n, B_1..B_n, C_1..C_n).
 *   Edge weights: A_i -> B_i = x_i, B_i -> C_i = x_i,
 *   C_i -> A_i = x_i (forming triangular gadgets).
 *   Then x_i + x_j + x_k = 0 iff there is a path of length
 *   exactly x_i + x_j + x_k = 0 from A_i to A_i through B_j
 *   and C_k, which can be detected via APSP.
 *
 * The reverse reduction (APSP -> 3SUM) encodes the min-plus
 * matrix product as a 3SUM-style problem over structured sets.
 *
 * For further reading:
 *   - Gajentaan & Overmars (1995): "On a class of O(n^2) problems
 *     in computational geometry", CGTA 5:165-185
 *   - Gronlund & Pettie (2014): "Threesomes, Degenerates, and
 *     Love Triangles", FOCS 2014
 *   - Gold & Sharir (2017): "Improved 3SUM Algorithm",
 *     STOC 2017
 *   - Vassilevska-Williams & Williams (2010): "Subcubic
 *     Equivalences between Path, Matrix, and Triangle Problems",
 *     FOCS 2010
 */
/*
 * Algorithm comparison:
 *
 * | Algorithm              | Time         | Space  | Deterministic |
 * |------------------------|--------------|--------|---------------|
 * | Naive O(n^3)           | Theta(n^3)   | O(1)   | Yes           |
 * | Hash O(n^2) expected   | O(n^2) avg   | O(n)   | No (random)   |
 * | Sort + BSearch         | O(n^2 log n) | O(1)   | Yes           |
 * | Two-Pointer (sorted)   | O(n^2)       | O(1)   | Yes           |
 * | Gronlund-Pettie (2014) | O(n^2/f(n))  | O(n)   | Yes           |
 * | Gold-Sharir (2017)     | O(n^2/g(n))  | O(n)   | Yes           |
 *
 * where f(n) = (log n / log log n)^{2/3} and g(n) = log n / log log n.
 *
 * All algorithms above achieve the same asymptotic n^2 behavior
 * up to polylogarithmic factors. The 3SUM Conjecture asserts that
 * no O(n^{2-epsilon}) algorithm exists for any epsilon > 0.
 */