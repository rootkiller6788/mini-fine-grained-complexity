/* condlb.h -- Conditional Lower Bounds: Core Data Structures and API
 *
 * Defines the formal framework for conditional lower bounds (CLB) in
 * fine-grained complexity theory. A conditional lower bound states:
 *   "If Hypothesis H holds, then Problem P requires time T(n)."
 *
 * This methodology replaces unconditional lower bounds (which are
 * notoriously hard to prove, even for P vs NP) with conditional
 * statements under plausible hardness assumptions.
 *
 * Key hypotheses in fine-grained complexity:
 *   SETH — Strong Exponential Time Hypothesis (Impagliazzo-Paturi, 2001)
 *   ETH  — Exponential Time Hypothesis (Impagliazzo-Paturi, 1999)
 *   OVC  — Orthogonal Vectors Conjecture (Williams, 2005)
 *   3SUM — 3SUM Conjecture (Gajentaan-Overmars, 1995; Patrascu, 2010)
 *   APSP — All-Pairs Shortest Paths Conjecture (Williams-Williams, 2010)
 *
 * References:
 *   Williams (2018): "Fine-Grained Complexity" (textbook)
 *   Williams (2015): "Hardness of Easy Problems: Basing Hardness on
 *                     Popular Conjectures such as SETH"
 *   Patrascu (2010): "Towards Polynomial Lower Bounds for Dynamic Problems"
 *   Abboud-Williams (2014): "Popular Conjectures Imply Strong Lower
 *                            Bounds for Dynamic Problems"
 */

#ifndef CONDLB_H
#define CONDLB_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <assert.h>

/* ============================================================================
 * L1: Core Definitions — Conditional Lower Bound Framework
 * ============================================================================ */

typedef enum {
    HYPOTHESIS_UNREFUTED  = 0,
    HYPOTHESIS_REFUTED    = 1,
    HYPOTHESIS_OPEN       = 2,
    HYPOTHESIS_EQUIVALENT = 3
} HypothesisStatus;

typedef struct Hypothesis {
    char             name[128];
    char             problem[128];
    double           lower_bound_exponent;
    double           epsilon_tolerance;
    HypothesisStatus status;
    int              year_proposed;
    char             proposed_by[128];
    int              n_depends_on;
    int*             depends_on_ids;
    int              n_implies;
    int*             implies_ids;
} Hypothesis;

typedef enum {
    REDUCTION_MANY_ONE     = 0,
    REDUCTION_TURING       = 1,
    REDUCTION_RANDOMIZED   = 2,
    REDUCTION_FINE_GRAINED = 3
} ReductionType;

typedef struct {
    int           hypothesis_id;
    char          target_problem[128];
    double        time_bound_exponent;
    ReductionType reduction_type;
    int           tightness;
    int           proven;
    char          reference[256];
} ConditionalLowerBound;

typedef struct {
    Hypothesis*            hypotheses;
    int                    n_hypotheses;
    int                    max_hypotheses;
    ConditionalLowerBound* bounds;
    int                    n_bounds;
    int                    max_bounds;
} HypothesisDatabase;

typedef struct {
    double exponent;
    double polylog_power;
    int    is_upper_bound;
    int    is_exponential;
    double exponential_base;
} FineGrainedTimeBound;

typedef struct {
    int    from_problem_id;
    int    to_problem_id;
    double time_blowup;
    double polylog_blowup;
    int    randomized;
    char   paper_ref[256];
} ReductionEdge;

typedef struct {
    char*            problem_names;
    int*             name_offsets;
    int              n_problems;
    int              max_problems;
    ReductionEdge*   edges;
    int              n_edges;
    int              max_edges;
    int**            adjacency;
    int*             adj_sizes;
    int*             adj_capacities;
} ReductionWeb;

/* API: Hypothesis Database Management */
HypothesisDatabase* hyp_db_create(int max_hypotheses, int max_bounds);
void                hyp_db_free(HypothesisDatabase* db);
int  hyp_db_add_hypothesis(HypothesisDatabase* db, const char* name,
                            const char* problem, double exponent,
                            double eps_tol, HypothesisStatus status,
                            int year, const char* proposed_by);
void hyp_db_add_implication(HypothesisDatabase* db, int from_id, int to_id);
int  hyp_db_add_bound(HypothesisDatabase* db, int hyp_id,
                      const char* target_problem, double exponent,
                      ReductionType rtype, int tightness, const char* reference);
int  hyp_db_find_hypothesis(const HypothesisDatabase* db, const char* name);
void hyp_db_describe(const HypothesisDatabase* db, int hyp_id, char* buf, int buf_sz);
void hyp_db_print(const HypothesisDatabase* db);

/* API: Reduction Web */
ReductionWeb* rw_create(int max_problems);
void          rw_free(ReductionWeb* rw);
int  rw_add_problem(ReductionWeb* rw, const char* name);
void rw_add_reduction(ReductionWeb* rw, int from_id, int to_id,
                      double time_blowup, double polylog_blowup,
                      int randomized, const char* paper_ref);
int  rw_path_exists(const ReductionWeb* rw, int from_id, int to_id);
int** rw_transitive_closure(const ReductionWeb* rw, int* n);
int  rw_equivalence_class(const ReductionWeb* rw, int problem_id,
                           int* output, int max_output);
int  rw_hardness_tier(const ReductionWeb* rw, int problem_id);
void rw_print(const ReductionWeb* rw);

/* API: Time Bounds */
FineGrainedTimeBound fgtb_make(double exponent, double polylog_power,
                                int is_upper, int is_exp, double exp_base);
int  fgtb_compare(FineGrainedTimeBound a, FineGrainedTimeBound b);
int  fgtb_violates(FineGrainedTimeBound b, ConditionalLowerBound* lb);
void fgtb_format(FineGrainedTimeBound tb, char* buf, int buf_sz);

#endif /* CONDLB_H */
