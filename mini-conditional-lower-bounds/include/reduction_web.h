/* reduction_web.h -- Reduction Web and Equivalence Classes
 *
 * The reduction web is the graph of known fine-grained reductions
 * between problems. Nodes are computational problems; directed edges
 * represent fine-grained reductions (A -> B means A reduces to B,
 * so hardness of A implies hardness of B).
 *
 * Equivalence classes group problems that are mutually reducible
 * under fine-grained reductions. These form the "subcubic" and
 * "subquadratic" equivalence classes in fine-grained complexity.
 *
 * Key equivalence classes:
 *   - Subcubic equivalence class: APSP, Negative Triangle,
 *     Radius, Median, Betweenness Centrality, etc.
 *   - Subquadratic equivalence class: OV, Edit Distance, LCS,
 *     Diameter, Frechet Distance, Dynamic Time Warping, etc.
 *   - 3SUM-hard class: 3SUM, 3-Points-on-a-Line, Polygon Containment
 *
 * References:
 *   Williams & Williams (2010): "Subcubic Equivalences Between
 *        Path, Matrix, and Triangle Problems"
 *   Abboud, Williams, Yu (2015): "More Applications of the
 *        Polynomial Method to Algorithm Design"
 */

#ifndef REDUCTION_WEB_H
#define REDUCTION_WEB_H

#include "condlb.h"
#include <stdint.h>

/* ============================================================================
 * L1: Definitions - Problem Nodes and Equivalence Classes
 * ============================================================================ */

/* ProblemCategory: high-level grouping of problems */
typedef enum {
    CAT_STRING       = 0, /* string algorithms (Edit, LCS, etc.) */
    CAT_GRAPH        = 1, /* graph algorithms (APSP, Diameter, etc.) */
    CAT_GEOMETRY     = 2, /* computational geometry (3SUM, etc.) */
    CAT_MATRIX       = 3, /* matrix problems (MM, etc.) */
    CAT_SATISFIABILITY = 4, /* SAT and CSP variants */
    CAT_DYNAMIC      = 5, /* dynamic data structure problems */
    CAT_OTHER        = 6
} ProblemCategory;

/* A problem node in the reduction web */
typedef struct {
    char            name[128];
    ProblemCategory category;
    double          naive_exponent;     /* exponent of naive algorithm */
    double          best_known_exponent;/* best known algorithm exponent */
    double          conditional_lower_bound; /* under strongest applicable hypothesis */
    int             is_hard;            /* 1 if CLB proven for this problem */
    int             is_complete;        /* 1 if complete for some class */
} ProblemNode;

/* EquivalenceClass: a set of mutually reducible problems */
typedef struct {
    int    class_id;
    char   representative[128];  /* name of the canonical problem */
    int*   member_ids;
    int    n_members;
    int    max_members;
    double hardness_exponent;    /* the exponent shared by all members */
    int    is_subcubic;          /* 1 if part of subcubic equivalence */
    int    is_subquadratic;      /* 1 if part of subquadratic equivalence */
} EquivalenceClass;

/* ============================================================================
 * L2: Core Concepts - Extended Reduction Web
 * ============================================================================ */

/* Extended reduction web with problem metadata and equivalence classes */
typedef struct {
    ProblemNode*      problems;
    int               n_problems;
    int               max_problems;
    ReductionEdge*    edges;
    int               n_edges;
    int               max_edges;
    int**             adjacency;
    int*              adj_sizes;
    int*              adj_capacities;
    EquivalenceClass* classes;
    int               n_classes;
    int               max_classes;
} ExtendedReductionWeb;

/* ============================================================================
 * L3: Mathematical Structures - Weighted Reductions
 * ============================================================================ */

/* A weighted reduction with exact time transformation.
 * If A reduces to B with time T_B(n) per oracle call and k calls,
 * then T_A(n) <= k * T_B(n^c) + overhead_O(n^d). */
typedef struct {
    int    from_id;
    int    to_id;
    double oracle_calls;      /* number of oracle calls k */
    double input_blowup;      /* transform input size n -> n^c */
    double overhead_exponent; /* d where overhead is O(n^d) */
} WeightedReduction;

/* LowerBoundChain: a sequence of reductions proving a CLB */
typedef struct {
    int              length;       /* number of reductions in the chain */
    WeightedReduction* reductions;
    double           final_exponent; /* the derived lower bound */
    int              valid;         /* 1 if chain is logically sound */
} LowerBoundChain;

/* ============================================================================
 * API: Extended Reduction Web
 * ============================================================================ */

ExtendedReductionWeb* erw_create(int max_problems, int max_edges, int max_classes);
void erw_free(ExtendedReductionWeb* erw);

/* Register a problem in the web. Returns problem id. */
int erw_add_problem(ExtendedReductionWeb* erw, const char* name,
                    ProblemCategory cat, double naive_exp,
                    double best_known_exp);

/* Add a reduction edge between problems. */
void erw_add_edge(ExtendedReductionWeb* erw, int from_id, int to_id,
                  double time_blowup, double polylog_blowup,
                  int randomized, const char* paper_ref);

/* Set the conditional lower bound for a problem. */
void erw_set_clb(ExtendedReductionWeb* erw, int problem_id,
                 double lower_bound);

/* Mark a problem as complete for its equivalence class. */
void erw_mark_complete(ExtendedReductionWeb* erw, int problem_id);

/* ============================================================================
 * API: Equivalence Class Computation
 * ============================================================================ */

/* Compute all equivalence classes (strongly connected components)
 * of the reduction web. Returns the number of classes found.
 * Uses Kosaraju's algorithm for SCC decomposition.
 * Time: O(V + E). */
int erw_compute_classes(ExtendedReductionWeb* erw);

/* Find the equivalence class containing a given problem.
 * Returns class id, or -1 if not found. */
int erw_find_class(const ExtendedReductionWeb* erw, int problem_id);

/* Get all problems in a given equivalence class.
 * Returns number of members written to output. */
int erw_class_members(const ExtendedReductionWeb* erw, int class_id,
                      int* output, int max_output);

/* ============================================================================
 * API: Lower Bound Chain Construction
 * ============================================================================ */

/* Construct a lower bound chain from a hypothesis to a target problem.
 * Uses BFS to find the shortest reduction path.
 * Returns a chain (caller must free with lbc_free).
 * Time: O(V + E). */
LowerBoundChain* lbc_construct(const ExtendedReductionWeb* erw,
                               int hypothesis_problem_id,
                               int target_problem_id);

/* Validate a lower bound chain by checking all reduction
 * parameters are consistent. Returns 1 if valid. */
int lbc_validate(const LowerBoundChain* chain);

/* Compute the final exponent implied by the chain.
 * Propagates time bounds through all reductions. */
double lbc_compute_exponent(const LowerBoundChain* chain);

void lbc_free(LowerBoundChain* chain);

/* ============================================================================
 * API: Subcubic and Subquadratic Equivalences
 * ============================================================================ */

/* Check if two problems are in the same subcubic equivalence class. */
int erw_same_subcubic_class(const ExtendedReductionWeb* erw,
                            int prob1_id, int prob2_id);

/* Check if two problems are in the same subquadratic equivalence class. */
int erw_same_subquadratic_class(const ExtendedReductionWeb* erw,
                                int prob1_id, int prob2_id);

/* Print the subcubic equivalence class with all members, reductions,
 * and the common hardness exponent. */
void erw_print_subcubic_class(const ExtendedReductionWeb* erw);

/* Print the subquadratic equivalence class. */
void erw_print_subquadratic_class(const ExtendedReductionWeb* erw);

/* Print a formatted summary of the entire reduction web. */
void erw_print_summary(const ExtendedReductionWeb* erw);

/* ============================================================================
 * API: Hardness Propagation
 * ============================================================================ */

/* Propagate hardness: given that a set of "source" hypotheses are
 * assumed true, compute the strongest conditional lower bound for
 * every reachable problem in the web.
 * Time: O(V * E) or O(V^3) for full propagation. */
void erw_propagate_hardness(ExtendedReductionWeb* erw);

/* Get the strongest lower bound that can be derived for a problem
 * given the current state of the reduction web. */
double erw_strongest_bound(const ExtendedReductionWeb* erw, int problem_id);

#endif /* REDUCTION_WEB_H */
