#ifndef OV_H
#define OV_H
#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* mini-orthogonal-vectors - OV Problem and SETH Connection
 * OV: Given A,B of n vectors in {0,1}^d, is there a.b=0?
 * OV Conjecture: n^{2-o(1)} for d=omega(log n). SETH=>OV (Williams 2005).
 */

/* L1: Binary vector in packed word representation (32 bits/word) */
typedef struct { uint32_t *bits; int32_t dimension, id; } binary_vector_t;
/* L1: Set of n vectors sharing dimension d */
typedef struct { binary_vector_t *vectors; int32_t num_vectors, dimension; } vector_set_t;
/* L1: OV instance = two vector sets A and B */
typedef struct { vector_set_t *A, *B; } ov_instance_t;
/* L1: OV query result with statistics */
typedef struct { bool found; int32_t a_index, b_index; int64_t ops; double time_s; } ov_result_t;
/* L1: Density statistics for a vector set */
typedef struct { double mean_density, min_density, max_density, stddev_density; int32_t num_vectors, num_zero_vectors, num_ones_vectors; } density_stats_t;

/* === L1: Binary Vector Operations === */
binary_vector_t *bv_create(int32_t dim);
binary_vector_t *bv_copy(const binary_vector_t *src);
void bv_copy_into(const binary_vector_t *src, binary_vector_t *dst);
void bv_destroy(binary_vector_t *v);
bool bv_get(const binary_vector_t *v, int32_t pos);
void bv_set(binary_vector_t *v, int32_t pos, bool val);
void bv_flip(binary_vector_t *v, int32_t pos);
void bv_clear(binary_vector_t *v);
void bv_fill_ones(binary_vector_t *v);
void bv_set_random(binary_vector_t *v, double density);
bool bv_set_from_string(binary_vector_t *v, const char *s);
int32_t bv_dot_product(const binary_vector_t *a, const binary_vector_t *b);
bool bv_are_orthogonal(const binary_vector_t *a, const binary_vector_t *b);
int32_t bv_hamming_weight(const binary_vector_t *v);
int32_t bv_hamming_distance(const binary_vector_t *a, const binary_vector_t *b);
bool bv_orthogonal_on_subset(const binary_vector_t *a, const binary_vector_t *b,
                              const int32_t *pos, int32_t n);
double bv_jaccard(const binary_vector_t *a, const binary_vector_t *b);
double bv_cosine(const binary_vector_t *a, const binary_vector_t *b);
void bv_print(const binary_vector_t *v);
void bv_fprint(FILE *fp, const binary_vector_t *v);

/* === L1: Vector Set Operations === */
vector_set_t *vs_create(int32_t n, int32_t dim);
vector_set_t *vs_copy(const vector_set_t *src);
void vs_destroy(vector_set_t *vs);
binary_vector_t *vs_get(vector_set_t *vs, int32_t i);
void vs_set_random(vector_set_t *vs, double density);
int32_t vs_set_from_strings(vector_set_t *vs, const char **pats, int32_t n);
density_stats_t vs_density_stats(const vector_set_t *vs);
void vs_print(const vector_set_t *vs);

/* === L1: OV Instance Operations === */
ov_instance_t *ov_create(int32_t n, int32_t d);
ov_instance_t *ov_copy(const ov_instance_t *inst);
void ov_destroy(ov_instance_t *inst);
int64_t ov_num_pairs(const ov_instance_t *inst);

/* === L2: OV Conjecture and Hardness === */
bool ov_conjecture_violated(int32_t n, int32_t d, double exp);
double ov_lower_bound_exponent(int32_t n, int32_t d);
bool ov_is_superlogarithmic(int32_t n, int32_t d);
int32_t ov_regime_classify(int32_t n, int32_t d);

/* === L3: Packed Bit Operations === */
int32_t bv_num_words(int32_t dim);
void bv_word_bit_index(int32_t pos, int32_t *wi, uint32_t *bm);
void bv_and(const binary_vector_t *a, const binary_vector_t *b, binary_vector_t *r);
void bv_or(const binary_vector_t *a, const binary_vector_t *b, binary_vector_t *r);
void bv_xor(const binary_vector_t *a, const binary_vector_t *b, binary_vector_t *r);
void bv_not(const binary_vector_t *v, binary_vector_t *r);
bool bv_is_zero(const binary_vector_t *v);
bool bv_is_ones(const binary_vector_t *v);
bool bv_equals(const binary_vector_t *a, const binary_vector_t *b);
int32_t bv_popcount32(uint32_t w);
bool bv_is_subset(const binary_vector_t *v, const binary_vector_t *w);

/* === L4: Fundamental Laws === */
double ov_to_ksat_threshold(double ov_exp, int32_t n, int32_t d);
bool ov_sparsification_applies(int32_t n, int32_t d);
int32_t ov_equivalent_ksat_param(double ratio);
bool ov_equivalent_to_seth(int32_t n, int32_t d);

/* === L5: Algorithms === */
ov_result_t ov_brute_force(const ov_instance_t *inst);
ov_result_t ov_brute_force_early(const ov_instance_t *inst);
ov_result_t ov_brute_force_sparse_first(const ov_instance_t *inst);
ov_result_t ov_williams(const ov_instance_t *inst);
ov_result_t ov_williams_configurable(const ov_instance_t *inst, int32_t dp, int32_t tr);
ov_result_t ov_meet_in_middle(const ov_instance_t *inst);
ov_result_t ov_four_russians(const ov_instance_t *inst);
ov_result_t ov_via_matrix_mult(const ov_instance_t *inst);
ov_result_t ov_polynomial_method(const ov_instance_t *inst);
ov_result_t ov_light_sparse(const ov_instance_t *inst);
ov_result_t ov_lsh(const ov_instance_t *inst, int32_t L, int32_t k);
int64_t ov_count_orthogonal_pairs(const ov_instance_t *inst);
bool ov_has_k_orthogonal_pairs(const ov_instance_t *inst, int64_t K);
int64_t ov_find_all_orthogonal_pairs(const ov_instance_t *inst,
                                      int32_t *ai, int32_t *bi, int64_t max);

/* === L6: Canonical Problems - Reductions === */
typedef struct { char *text, *pattern; int32_t text_len, pattern_len; } pm_instance_t;
typedef struct { int32_t num_vertices, num_edges; int32_t *edges_u, *edges_v; int32_t *degree, *adj_offset, *adj_list; } sparse_graph_t;
typedef struct { char *str1, *str2; int32_t len1, len2; } ed_instance_t;
typedef struct { char *seq1, *seq2; int32_t len1, len2; } lcs_instance_t;
typedef struct { int64_t *numbers; int32_t count; int64_t target; } subset_sum_instance_t;

pm_instance_t *ov_reduce_to_pattern_matching(const ov_instance_t *inst);
void pm_instance_destroy(pm_instance_t *pm);
int32_t pm_find_naive(const pm_instance_t *pm, int32_t *m, int32_t maxm);
sparse_graph_t *ov_reduce_to_graph_diameter(const ov_instance_t *inst);
void sparse_graph_destroy(sparse_graph_t *g);
int32_t graph_diameter_naive(const sparse_graph_t *g);
ed_instance_t *ov_reduce_to_edit_distance(const ov_instance_t *inst);
void ed_instance_destroy(ed_instance_t *ed);
int32_t edit_distance_naive(const char *s1, int32_t l1, const char *s2, int32_t l2);
lcs_instance_t *ov_reduce_to_lcs(const ov_instance_t *inst);
void lcs_instance_destroy(lcs_instance_t *lcs);
int32_t lcs_length_naive(const char *s1, int32_t l1, const char *s2, int32_t l2);
subset_sum_instance_t *ov_reduce_to_subset_sum(const ov_instance_t *inst);
void subset_sum_instance_destroy(subset_sum_instance_t *ss);

/* === L7: Applications === */
int64_t ov_database_anti_join(const ov_instance_t *inst, int32_t *ra, int32_t *rb, int64_t max);
ov_instance_t *ov_from_dna_sequences(const char **seqs, int32_t n, int32_t len);
int32_t ov_count_safe_crispr_pairs(const ov_instance_t *inst);
double ov_max_code_dissimilarity(const ov_instance_t *inst, int32_t *fa, int32_t *fb);
int32_t ov_find_topk_diverse_users(const ov_instance_t *inst, int32_t k,
                                    int32_t *ua, int32_t *ub, double *sc);

/* === L8: Advanced Topics === */
int32_t ov_polynomial_degree(int32_t d);
int32_t *ov_polynomial_coefficients(int32_t d, int32_t *nc);
int64_t ov_communication_lower_bound(int32_t n, int32_t d);
int64_t ov_communication_protocol_simulate(const ov_instance_t *inst, int32_t bpm);
ov_result_t ov_streaming_online(const ov_instance_t *inst);

typedef struct {
    vector_set_t *A, *B; int32_t max_a, max_b, count_a, count_b;
    bool has_answer; int32_t cached_a, cached_b;
} ov_dynamic_t;

ov_dynamic_t *ov_dynamic_create(int32_t maxn, int32_t d);
bool ov_dynamic_insert_a(ov_dynamic_t *dyn, const binary_vector_t *v);
bool ov_dynamic_insert_b(ov_dynamic_t *dyn, const binary_vector_t *v);
bool ov_dynamic_query(ov_dynamic_t *dyn, int32_t *ai, int32_t *bi);
bool ov_dynamic_delete_a(ov_dynamic_t *dyn, int32_t i);
bool ov_dynamic_delete_b(ov_dynamic_t *dyn, int32_t i);
void ov_dynamic_destroy(ov_dynamic_t *dyn);
ov_result_t ov_approximate_orthogonal(const ov_instance_t *inst, double eps);

/* === L9: Research Frontiers === */
bool ov_conjecture_implies_3sum(int32_t n, int32_t d);
void ov_print_equivalence_class(void);

/* === Instance Generation === */
ov_instance_t *ov_create_random(int32_t n, int32_t d, double density, bool guarantee);
ov_instance_t *ov_create_from_ksat(int32_t nv, int32_t k, int32_t nc);
ov_instance_t *ov_create_worst_case(int32_t n, int32_t d);
ov_instance_t *ov_create_all_orthogonal(int32_t n, int32_t d);
ov_instance_t *ov_create_no_orthogonal(int32_t n, int32_t d);
ov_instance_t *ov_create_unique_pair(int32_t n, int32_t d, double density);

/* === Utilities === */
double ov_get_time_seconds(void);
void ov_print_summary(const ov_instance_t *inst);
bool ov_validate(const ov_instance_t *inst);
void ov_report_result(const ov_result_t *r);

/* === Utility: Benchmarking & Cross-Validation === */
bool ov_compare_algorithms(const ov_instance_t *inst,
                            const char *n1, ov_result_t (*a1)(const ov_instance_t*),
                            const char *n2, ov_result_t (*a2)(const ov_instance_t*));
ov_result_t ov_benchmark(const ov_instance_t *inst,
                          ov_result_t (*alg)(const ov_instance_t*), int32_t repeats);
void ov_benchmark_all(const ov_instance_t *inst);
double ov_estimate_exponent(int32_t d, double density, int32_t mn, int32_t mx, int32_t st);
double ov_orthogonal_fraction(const ov_instance_t *inst);
bool ov_verify_conjecture_empirical(int32_t n, int32_t d);
void ov_print_complexity_table(void);
int32_t *ov_bit_frequencies(const ov_instance_t *inst, bool set_a);
void ov_print_bit_histogram(const ov_instance_t *inst);
int32_t ov_self_test(void);

/* === Utility: Extended Operations === */
void bv_add_gf2(const binary_vector_t *a, const binary_vector_t *b, binary_vector_t *r);
void bv_scalar_mul_gf2(const binary_vector_t *v, bool s, binary_vector_t *r);
void bv_linear_combination_gf2(const binary_vector_t *a, bool c1,
                                const binary_vector_t *b, bool c2, binary_vector_t *r);
bool bv_linearly_independent_gf2(const binary_vector_t *a, const binary_vector_t *b);
int32_t vs_rank_gf2(const vector_set_t *vs);
int32_t ov_min_distance_to_orthogonal(const ov_instance_t *inst);
int32_t ov_max_overlap(const ov_instance_t *inst);
void ov_print_implication_chain(void);
double ov_instance_difficulty(const ov_instance_t *inst);
void ov_print_difficulty_report(const ov_instance_t *inst);
double ov_expected_orthogonal_fraction(int32_t n, int32_t d, double density);
double ov_expected_orthogonal_count(int32_t n, int32_t d, double density);
void ov_print_expected_vs_actual(const ov_instance_t *inst);
void ov_empirical_seth_check(int32_t d, double density, int32_t ns, int32_t ne, int32_t step);
void ov_print_ascii_visualization(const ov_instance_t *inst);

#endif /* OV_H */
