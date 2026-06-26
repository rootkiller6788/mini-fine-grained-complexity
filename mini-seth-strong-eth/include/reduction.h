#ifndef FINE_GRAINED_REDUCTION_H
#define FINE_GRAINED_REDUCTION_H

#include "seth.h"

/* ============================================================================
 * Fine-Grained Reductions
 *
 * L2: Core concepts of subcubic and subquadratic reductions
 * L4: Fine-grained reduction completeness
 *
 * A fine-grained reduction from problem A to problem B shows:
 * If B can be solved in time T(n), then A can be solved in time T(n)^{O(1)}.
 * The key is preserving the exact exponent.
 * ============================================================================ */

/* ---- Reduction Types ---- */

typedef enum {
    REDUCTION_KARP,            /* standard Karp (many-one) reduction */
    REDUCTION_TURING,          /* Turing (oracle) reduction */
    REDUCTION_SUBCUBIC,        /* preserves subcubic equivalence */
    REDUCTION_SUBQUADRATIC,    /* preserves subquadratic equivalence */
    REDUCTION_FINE_GRAINED     /* general fine-grained reduction */
} reduction_type_t;

/* A fine-grained reduction record */
typedef struct {
    reduction_type_t type;
    char            *from_problem;
    char            *to_problem;
    double           time_factor;   /* multiplicative overhead */
    double           size_blowup;   /* n → n' = O(n^c) */
    char            *description;
} fg_reduction_t;

/* ---- Reduction Graph ---- */

/* Directed graph of fine-grained reductions between problems */
typedef struct {
    char          **problem_names;
    int32_t         num_problems;
    int32_t         capacity;
    fg_reduction_t **adjacency;  /* adjacency[i][j] = reduction i→j or NULL */
} reduction_graph_t;

reduction_graph_t *reduction_graph_create(int32_t initial_capacity);
void reduction_graph_destroy(reduction_graph_t *g);
int32_t reduction_graph_add_problem(reduction_graph_t *g, const char *name);
void reduction_graph_add_edge(reduction_graph_t *g,
                               int32_t from, int32_t to,
                               const fg_reduction_t *r);
fg_reduction_t *reduction_graph_get_edge(const reduction_graph_t *g,
                                          int32_t from, int32_t to);

/* ---- Specific Fine-Grained Reductions ---- */

/* SETH → OV reduction:
 * If OV_{n,d} can be solved in O(n^{2-ε}) time for d=ω(log n),
 * then SETH is false.
 * Reference: Williams (2005) */
typedef struct {
    int32_t num_vectors;
    int32_t dimension;
    int32_t *bits;     /* packed bits: num_vectors * dimension / 32 + 1 */
} ov_instance_t;

/* Reduce k-SAT to OV instance */
ov_instance_t *seth_to_ov_reduction(const cnf_formula_t *f, int32_t k);

/* SETH → Edit Distance reduction:
 * If Edit Distance can be solved in O(n^{2-ε}) time, SETH is false.
 * Reference: Backurs, Indyk (2015) */
typedef struct {
    char   *s1;
    char   *s2;
    int32_t n;
} edit_distance_instance_t;

edit_distance_instance_t *seth_to_edit_distance(const cnf_formula_t *f);

/* SETH → k-Clique reduction:
 * If k-Clique (small k) can be solved faster than exhaustive search, SETH false.
 * Reference: Abboud, Williams, Yu (2015) */

/* SETH → LCS reduction:
 * If LCS can be solved in O(n^{2-ε}) time, SETH is false.
 * Reference: Abboud, Backurs, Williams (2015) */

/* ---- Reduction Chain Utilities ---- */

/* Check if a reduction chain from SETH to problem P exists */
bool reduction_chain_exists(const reduction_graph_t *g,
                             const char *from, const char *to);

/* Compute the conditional lower bound for a problem based on SETH */
double seth_lower_bound(const reduction_graph_t *g, const char *problem);

/* Print the reduction path */
void print_reduction_path(const reduction_graph_t *g,
                           const char *problem);

/* Build the standard reduction graph with all known reductions */
reduction_graph_t *build_standard_reduction_graph(void);

#endif /* FINE_GRAINED_REDUCTION_H */
