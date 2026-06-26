#ifndef EQUIV_CLASSES_H
#define EQUIV_CLASSES_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <math.h>

/*
 * mini-equivalence-classes: Fine-Grained Equivalence Classes of Problems
 *
 * Core definitions for fine-grained equivalence classes:
 *   - Subcubic equivalence class (APSP, Negative Triangle, Radius, Median)
 *   - Subquadratic equivalence class (OV, Edit Distance, LCS, Frechet, DTW)
 *   - 3SUM equivalence class
 *   - Reduction-preserving equivalences
 *
 * References:
 *   Williams (2015), "Hardness of Easy Problems", ICDT
 *   Williams & Williams (2013), STOC 2013
 *   Backurs & Indyk (2016), SIAM J. Computing
 *   Gajentaan & Overmars (1995), Computational Geometry
 *
 * L1: Core type definitions for equivalence classes and problem instances
 * L2: Fine-grained reduction concepts, equivalence relation properties
 * L3: Mathematical structures: distance matrices, vector sets, point sets
 */

typedef enum {
    EQUIV_CLASS_SUBCUBIC       = 0,
    EQUIV_CLASS_SUBQUADRATIC   = 1,
    EQUIV_CLASS_3SUM           = 2,
    EQUIV_CLASS_APSP           = 3,
    EQUIV_CLASS_CNFSAT         = 4,
    EQUIV_CLASS_ALLPAIRS_LCS   = 5,
    EQUIV_CLASS_BMM            = 6,
    EQUIV_CLASS_DIAMETER       = 7,
    EQUIV_CLASS_HITTING_SET    = 8,
    EQUIV_CLASS_COLLINEARITY   = 9,
    EQUIV_CLASS_COUNT          = 10
} equiv_class_id_t;

const char *equiv_class_name(equiv_class_id_t id);

typedef enum {
    FGR_SUBCUBIC     = 0,
    FGR_SUBQUADRATIC = 1,
    FGR_3SUM         = 2,
    FGR_LINEAR       = 3,
    FGR_NEAR_LINEAR  = 4,
    FGR_MANY_ONE     = 5,
    FGR_TURING       = 6
} fg_reduction_type_t;

typedef uint32_t problem_id_t;

typedef struct {
    problem_id_t        problem_id;
    equiv_class_id_t    equiv_class;
    char                name[80];
    int32_t             input_size;
    int32_t             secondary_size;
    double              time_exponent;
    bool                is_conjectured_hard;
} problem_descriptor_t;

typedef struct {
    problem_id_t        from_problem;
    problem_id_t        to_problem;
    fg_reduction_type_t reduction_type;
    double              runtime_exponent;
    int32_t             output_blowup;
    bool                is_bi_reduction;
    bool                preserves_exact_exponent;
    char                citation[160];
} reduction_descriptor_t;

typedef struct {
    equiv_class_id_t    class_id;
    char                canonical_name[64];
    double              threshold_exponent;
    int32_t             num_members;
    int32_t             num_reductions;
    double              best_algorithm_exponent;
    double              conjectured_lower_bound;
    bool                is_runtime_equivalent;
} equiv_class_t;

equiv_class_t *equiv_class_create(equiv_class_id_t id, const char *canonical_name,
                                   double threshold, double conjectured_lb);
void equiv_class_free(equiv_class_t *ec);
bool equiv_class_add_member(equiv_class_t *ec, problem_id_t problem_id);
bool equiv_class_contains(const equiv_class_t *ec, problem_id_t problem_id);
bool problems_in_same_class(problem_id_t p1, problem_id_t p2);
bool register_reduction(problem_id_t from, problem_id_t to,
                        fg_reduction_type_t type, double exponent,
                        int32_t blowup, const char *citation);
double reduction_chain_exponent(problem_id_t from, problem_id_t to);
double class_diameter_exponent(equiv_class_id_t class_id);
bool class_is_tight(equiv_class_id_t class_id);
bool is_canonical_problem(problem_id_t problem_id);
problem_id_t class_canonical_problem(equiv_class_id_t class_id);
bool would_refute_subcubic_conjecture(double runtime_exponent, int32_t n);
bool would_refute_ov_conjecture(double runtime_exponent, int32_t n, int32_t d);
bool would_refute_3sum_conjecture(double runtime_exponent, int32_t n);
double effective_exponent_subcubic(int32_t n, double omega);
equiv_class_id_t classify_problem_by_structure(const problem_descriptor_t *desc);
int compare_equivalence_classes(equiv_class_id_t a, equiv_class_id_t b);
void print_equivalence_class_summary(void);
void print_class_reduction_graph(equiv_class_id_t class_id);
bool class_a_at_least_as_hard_as_b(equiv_class_id_t a, equiv_class_id_t b,
                                    const char *assumption);
double fg_landscape_distance(equiv_class_id_t a, equiv_class_id_t b);
bool class_a_implies_class_b(equiv_class_id_t a, equiv_class_id_t b);
problem_id_t register_problem(const char *name, equiv_class_id_t default_class,
                               double best_exponent, int32_t size_param);
const problem_descriptor_t *get_problem_descriptor(problem_id_t id);
int32_t total_registered_problems(void);
void print_problem_info(problem_id_t problem_id);

#endif /* EQUIV_CLASSES_H */