#ifndef OV_REDUCTIONS_H
#define OV_REDUCTIONS_H
#include "ov.h"

/* Reduction function declarations for L6 canonical problems.
 * All types (pm_instance_t, sparse_graph_t, etc.) are defined in ov.h. */

pm_instance_t *ov_reduce_to_pattern_matching(const ov_instance_t *inst);
void pm_instance_destroy(pm_instance_t *pm);
int32_t pm_find_naive(const pm_instance_t *pm, int32_t *matches, int32_t max_matches);
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
#endif
