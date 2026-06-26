#ifndef EDIT_DISTANCE_H
#define EDIT_DISTANCE_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <math.h>

/* ============================================================================
 * mini-edit-distance-lcs-hardness: Edit Distance
 *
 * Core definitions and data structures for edit distance (Levenshtein
 * distance) computation and conditional lower bounds.
 *
 * Reference: Levenshtein (1966) "Binary codes capable of correcting
 *            deletions, insertions, and reversals."
 *            Wagner & Fischer (1974) "The String-to-String Correction Problem"
 *            Backurs & Indyk (STOC 2015) "Edit Distance Cannot Be Computed
 *            in Strongly Subquadratic Time (unless SETH is false)"
 *
 * L1: Core definitions: edit distance, operations, DP table
 * L2: Quadratic-time lower bound under SETH
 * L3: DP matrix structure and operations
 * ============================================================================ */

/* ---- L1: Core Type Definitions ---- */

/* Edit operation types */
typedef enum {
    EDIT_OP_MATCH     = 0,
    EDIT_OP_SUBSTITUTE = 1,
    EDIT_OP_INSERT    = 2,
    EDIT_OP_DELETE    = 3
} edit_op_type_t;

/* A single edit operation */
typedef struct {
    edit_op_type_t type;
    int32_t        source_pos;
    int32_t        target_pos;
    char           source_char;
    char           target_char;
    int32_t        cost;
} edit_op_t;

/* An edit script: sequence of operations to transform source to target */
typedef struct {
    edit_op_t  *operations;
    int32_t     num_operations;
    int32_t     total_cost;
} edit_script_t;

/* Cost model for edit operations */
typedef struct {
    int32_t insert_cost;
    int32_t delete_cost;
    int32_t substitute_cost;
    int32_t match_cost;
} edit_cost_model_t;

/* ---- L2: Edit Distance Computation ---- */

int32_t edit_distance_full(const char *src, const char *tgt);
int32_t edit_distance_weighted(const char *src, const char *tgt,
                               const edit_cost_model_t *costs);
int32_t edit_distance_linear_space(const char *src, const char *tgt);
int32_t edit_distance_bounded(const char *src, const char *tgt, int32_t k);

/* ---- L3: DP Table Operations ---- */

typedef struct {
    int32_t  *data;
    int32_t   rows;
    int32_t   cols;
    int32_t   capacity;
} edit_dp_table_t;

typedef struct {
    int32_t       value;
    edit_op_type_t backpointer;
} edit_dp_cell_t;

typedef struct {
    edit_dp_cell_t *data;
    int32_t         rows;
    int32_t         cols;
} edit_dp_traceable_t;

edit_dp_table_t *ed_dp_create(int32_t src_len, int32_t tgt_len);
void ed_dp_destroy(edit_dp_table_t *dp);
void ed_dp_fill(edit_dp_table_t *dp, const char *src, const char *tgt);

edit_dp_traceable_t *ed_dp_traceable_create(int32_t src_len, int32_t tgt_len);
void ed_dp_traceable_destroy(edit_dp_traceable_t *dp);
void ed_dp_traceable_fill(edit_dp_traceable_t *dp, const char *src,
                           const char *tgt, const edit_cost_model_t *costs);

/* ---- Edit Script Operations ---- */

edit_script_t *ed_reconstruct_script(const edit_dp_traceable_t *dp,
                                     const char *src, const char *tgt);
edit_script_t *ed_reconstruct_linear(const char *src, const char *tgt);
void ed_script_destroy(edit_script_t *script);
void ed_script_print(const edit_script_t *script);

/* ---- Hamming Distance ---- */

int32_t hamming_distance(const char *s1, const char *s2);
double hamming_distance_weighted(const char *s1, const char *s2,
                                 const double *weights);

/* ---- Distance Metrics ---- */

int32_t damerau_levenshtein(const char *src, const char *tgt);
double jaro_similarity(const char *s1, const char *s2);
double jaro_winkler(const char *s1, const char *s2, double prefix_scale);
int32_t longest_common_prefix(const char *s1, const char *s2);
int32_t longest_common_suffix(const char *s1, const char *s2);
double edit_distance_normalized(const char *src, const char *tgt);
edit_script_t *edit_distance_with_script(const char *src, const char *tgt);

/* ---- L5: Faster Edit Distance Algorithms ---- */

int32_t edit_distance_four_russians(const char *src, const char *tgt);
int32_t edit_distance_banded(const char *src, const char *tgt, int32_t w);

/* ---- L8: Approximate Edit Distance ---- */

double edit_distance_approximate(const char *src, const char *tgt,
                                  double epsilon);
bool edit_distance_below_threshold_estimate(const char *src, const char *tgt,
                                             int32_t threshold);

#endif /* EDIT_DISTANCE_H */
