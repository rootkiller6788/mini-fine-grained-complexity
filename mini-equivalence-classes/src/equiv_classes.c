/* ============================================================================
 * equiv_classes.c -- Core Equivalence Class Framework Implementation
 *
 * Implements the fundamental data structures and operations for working
 * with fine-grained equivalence classes of computational problems.
 *
 * Knowledge Coverage:
 *   L1: Equivalence class creation/destruction, problem registration
 *   L2: Fine-grained reduction registration, equivalence checking
 *   L3: Internal data structures (hash tables, arrays, linked lists)
 *   L4: Reflexivity, transitivity, symmetry of equivalence classes
 *   L5: Reduction chain computation (Dijkstra on reduction graph)
 *   L7: Problem classification, comparison of equivalence classes
 *   L8: Cross-class relations, fine-grained landscape analysis
 *
 * References:
 *   Williams (2015), "Hardness of Easy Problems", ICDT
 *   Williams & Williams (2013), STOC 2013
 * ============================================================================ */

#include "equiv_classes.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ---- Internal Constants ---- */

#define MAX_PROBLEMS      256
#define MAX_REDUCTIONS    1024
#define MAX_CLASSES       16

/* ---- Internal Globals ---- */

static struct {
    problem_descriptor_t  problems[MAX_PROBLEMS];
    int32_t               num_problems;

    reduction_descriptor_t reductions[MAX_REDUCTIONS];
    int32_t                num_reductions;

    equiv_class_t         classes[MAX_CLASSES];
    bool                  class_initialized[MAX_CLASSES];

    /* Reduction graph adjacency matrix for path finding */
    double                adj[MAX_PROBLEMS][MAX_PROBLEMS];
    bool                  adj_valid[MAX_PROBLEMS][MAX_PROBLEMS];
} ctx;

static bool ctx_initialized = false;

/* ---- L1: Initialization ---- */

static void ctx_init(void) {
    if (ctx_initialized) return;
    memset(&ctx, 0, sizeof(ctx));
    for (int i = 0; i < MAX_PROBLEMS; i++) {
        for (int j = 0; j < MAX_PROBLEMS; j++) {
            ctx.adj[i][j] = -1.0;
            ctx.adj_valid[i][j] = false;
        }
    }
    for (int i = 0; i < MAX_CLASSES; i++) {
        ctx.class_initialized[i] = false;
    }
    ctx_initialized = true;
}

/* ---- L1: Equivalence Class Names ---- */

const char *equiv_class_name(equiv_class_id_t id) {
    static const char *names[] = {
        "Subcubic (APSP-equivalent)",
        "Subquadratic (OV-equivalent)",
        "3SUM-equivalent",
        "APSP",
        "CNF-SAT",
        "All-Pairs LCS",
        "Boolean Matrix Multiplication",
        "Graph Diameter/Radius",
        "Hitting Set / Set Cover",
        "Collinearity / 3SUM-hard",
        "UNKNOWN"
    };
    if (id < 0 || id >= EQUIV_CLASS_COUNT) return names[EQUIV_CLASS_COUNT];
    return names[id];
}

/* ---- L2: Problem Registration ---- */

problem_id_t register_problem(const char *name, equiv_class_id_t default_class,
                               double best_exponent, int32_t size_param) {
    ctx_init();
    if (ctx.num_problems >= MAX_PROBLEMS) return (problem_id_t)-1;
    problem_id_t id = (problem_id_t)ctx.num_problems;
    problem_descriptor_t *p = &ctx.problems[id];
    p->problem_id = id;
    p->equiv_class = default_class;
    strncpy(p->name, name, sizeof(p->name) - 1);
    p->name[sizeof(p->name) - 1] = '\0';
    p->input_size = size_param;
    p->secondary_size = 0;
    p->time_exponent = best_exponent;
    p->is_conjectured_hard = false;
    ctx.num_problems++;
    return id;
}

const problem_descriptor_t *get_problem_descriptor(problem_id_t id) {
    ctx_init();
    if (id >= ctx.num_problems) return NULL;
    return &ctx.problems[id];
}

int32_t total_registered_problems(void) {
    ctx_init();
    return ctx.num_problems;
}

void print_problem_info(problem_id_t problem_id) {
    ctx_init();
    if (problem_id >= ctx.num_problems) {
        printf("  Problem ID %u not found.\n", (unsigned)problem_id);
        return;
    }
    const problem_descriptor_t *p = &ctx.problems[problem_id];
    printf("  Problem [%u]: %s\n", (unsigned)problem_id, p->name);
    printf("    Class: %s\n", equiv_class_name(p->equiv_class));
    printf("    Best exponent: %.4f\n", p->time_exponent);
    printf("    Conjectured hard: %s\n", p->is_conjectured_hard ? "yes" : "no");
}

/* ---- L2: Equivalence Class Operations ---- */

equiv_class_t *equiv_class_create(equiv_class_id_t id, const char *canonical_name,
                                   double threshold, double conjectured_lb) {
    ctx_init();
    if (id < 0 || id >= MAX_CLASSES) return NULL;
    equiv_class_t *ec = &ctx.classes[id];
    memset(ec, 0, sizeof(equiv_class_t));
    ec->class_id = id;
    strncpy(ec->canonical_name, canonical_name, sizeof(ec->canonical_name) - 1);
    ec->canonical_name[sizeof(ec->canonical_name) - 1] = '\0';
    ec->threshold_exponent = threshold;
    ec->conjectured_lower_bound = conjectured_lb;
    ec->best_algorithm_exponent = threshold;
    ec->num_members = 0;
    ec->num_reductions = 0;
    ec->is_runtime_equivalent = false;
    ctx.class_initialized[id] = true;
    return ec;
}

void equiv_class_free(equiv_class_t *ec) {
    if (!ec) return;
    ctx.class_initialized[ec->class_id] = false;
    memset(ec, 0, sizeof(equiv_class_t));
}

bool equiv_class_add_member(equiv_class_t *ec, problem_id_t problem_id) {
    if (!ec || problem_id >= ctx.num_problems) return false;
    problem_descriptor_t *p = &ctx.problems[problem_id];
    p->equiv_class = ec->class_id;
    ec->num_members++;
    return true;
}

bool equiv_class_contains(const equiv_class_t *ec, problem_id_t problem_id) {
    if (!ec || problem_id >= ctx.num_problems) return false;
    return ctx.problems[problem_id].equiv_class == ec->class_id;
}

bool problems_in_same_class(problem_id_t p1, problem_id_t p2) {
    if (p1 >= ctx.num_problems || p2 >= ctx.num_problems) return false;
    return ctx.problems[p1].equiv_class == ctx.problems[p2].equiv_class;
}

/* ---- L2: Reduction Registration ---- */

bool register_reduction(problem_id_t from, problem_id_t to,
                        fg_reduction_type_t type, double exponent,
                        int32_t blowup, const char *citation) {
    ctx_init();
    if (from >= ctx.num_problems || to >= ctx.num_problems) return false;
    if (ctx.num_reductions >= MAX_REDUCTIONS) return false;

    reduction_descriptor_t *r = &ctx.reductions[ctx.num_reductions];
    r->from_problem = from;
    r->to_problem = to;
    r->reduction_type = type;
    r->runtime_exponent = exponent;
    r->output_blowup = blowup;
    r->is_bi_reduction = false;
    r->preserves_exact_exponent = (exponent < ctx.problems[from].time_exponent);
    strncpy(r->citation, citation, sizeof(r->citation) - 1);
    r->citation[sizeof(r->citation) - 1] = '\0';

    /* Update adjacency matrix */
    ctx.adj[from][to] = exponent;
    ctx.adj_valid[from][to] = true;

    /* Update class reduction count */
    equiv_class_id_t cid = ctx.problems[from].equiv_class;
    if (cid < MAX_CLASSES && ctx.class_initialized[cid]) {
        ctx.classes[cid].num_reductions++;
    }

    ctx.num_reductions++;
    return true;
}

/* ---- L2: Reduction Chain Computation (Floyd-Warshall) ---- */

double reduction_chain_exponent(problem_id_t from, problem_id_t to) {
    ctx_init();
    int32_t n = ctx.num_problems;
    if (from >= n || to >= n) return -1.0;
    if (from == to) return 0.0;

    /* Copy adjacency matrix */
    double dist[MAX_PROBLEMS][MAX_PROBLEMS];
    for (int32_t i = 0; i < n; i++) {
        for (int32_t j = 0; j < n; j++) {
            dist[i][j] = ctx.adj_valid[i][j] ? ctx.adj[i][j] : -1.0;
        }
        dist[i][i] = 0.0;
    }

    /* Floyd-Warshall for shortest path (min-plus semiring) */
    for (int32_t k = 0; k < n; k++) {
        for (int32_t i = 0; i < n; i++) {
            if (dist[i][k] < 0) continue;
            for (int32_t j = 0; j < n; j++) {
                if (dist[k][j] < 0) continue;
                double candidate = dist[i][k] + dist[k][j];
                if (dist[i][j] < 0 || candidate < dist[i][j]) {
                    dist[i][j] = candidate;
                }
            }
        }
    }

    return dist[from][to];
}

/* ---- L2: Equivalence Class Analysis ---- */

double class_diameter_exponent(equiv_class_id_t class_id) {
    ctx_init();
    double max_dist = 0.0;
    bool found_any = false;
    int32_t n = ctx.num_problems;

    for (int32_t i = 0; i < n; i++) {
        if (ctx.problems[i].equiv_class != class_id) continue;
        for (int32_t j = 0; j < n; j++) {
            if (ctx.problems[j].equiv_class != class_id) continue;
            double d = reduction_chain_exponent(i, j);
            if (d >= 0) {
                if (!found_any || d > max_dist) max_dist = d;
                found_any = true;
            }
        }
    }
    return found_any ? max_dist : -1.0;
}

bool class_is_tight(equiv_class_id_t class_id) {
    ctx_init();
    int32_t n = ctx.num_problems;
    double threshold = 0.0;
    bool found = false;

    /* Find threshold exponent for this class */
    for (int32_t i = 0; i < n; i++) {
        if (ctx.problems[i].equiv_class == class_id) {
            threshold = ctx.problems[i].time_exponent;
            found = true;
            break;
        }
    }
    if (!found) return false;

    /* Check: every pair has bidirectional reduction below threshold */
    for (int32_t i = 0; i < n; i++) {
        if (ctx.problems[i].equiv_class != class_id) continue;
        for (int32_t j = i + 1; j < n; j++) {
            if (ctx.problems[j].equiv_class != class_id) continue;
            double dij = reduction_chain_exponent(i, j);
            double dji = reduction_chain_exponent(j, i);
            if (dij < 0 || dij >= threshold) return false;
            if (dji < 0 || dji >= threshold) return false;
        }
    }
    return true;
}

bool is_canonical_problem(problem_id_t problem_id) {
    ctx_init();
    if (problem_id >= ctx.num_problems) return false;
    equiv_class_id_t cid = ctx.problems[problem_id].equiv_class;
    if (!ctx.class_initialized[cid]) return false;
    /* Canonical = name matches the class canonical name */
    const char *canonical = ctx.classes[cid].canonical_name;
    return strcmp(ctx.problems[problem_id].name, canonical) == 0;
}

problem_id_t class_canonical_problem(equiv_class_id_t class_id) {
    ctx_init();
    int32_t n = ctx.num_problems;
    const char *canonical = ctx.classes[class_id].canonical_name;
    for (int32_t i = 0; i < n; i++) {
        if (ctx.problems[i].equiv_class == class_id &&
            strcmp(ctx.problems[i].name, canonical) == 0) {
            return (problem_id_t)i;
        }
    }
    return (problem_id_t)-1;
}

/* ---- L2: Conjecture Verification ---- */

bool would_refute_subcubic_conjecture(double runtime_exponent, int32_t n) {
    (void)n;
    /* Subcubic conjecture: APSP cannot be solved in O(n^{3-eps}) */
    return runtime_exponent < 2.999;
}

bool would_refute_ov_conjecture(double runtime_exponent, int32_t n, int32_t d) {
    /* OV conjecture: no O(n^{2-eps}) when d = omega(log n) */
    if (d < (int32_t)(log2((double)n) * 2.0)) return false;
    return runtime_exponent < 1.999;
}

bool would_refute_3sum_conjecture(double runtime_exponent, int32_t n) {
    (void)n;
    /* 3SUM conjecture: no O(n^{2-eps}) for 3SUM */
    return runtime_exponent < 1.999;
}

double effective_exponent_subcubic(int32_t n, double omega) {
    /* Compute the effective exponent for subcubic problems */
    double log_n = log2((double)n);
    double mm_time = pow((double)n, omega); /* Fast MM */
    double mm_exponent = log2(mm_time) / log_n;
    /* Best of Floyd-Warshall and fast MM */
    double fw_exponent = 3.0;
    double best = (mm_exponent < fw_exponent) ? mm_exponent : fw_exponent;
    return best;
}

/* ---- L7: Problem Classification ---- */

equiv_class_id_t classify_problem_by_structure(const problem_descriptor_t *desc) {
    if (!desc) return EQUIV_CLASS_COUNT;
    double e = desc->time_exponent;
    if (e >= 2.9 && e <= 3.1) return EQUIV_CLASS_SUBCUBIC;
    if (e >= 1.9 && e <= 2.1) {
        if (strstr(desc->name, "3SUM") || strstr(desc->name, "Collinear"))
            return EQUIV_CLASS_3SUM;
        return EQUIV_CLASS_SUBQUADRATIC;
    }
    if (e >= 1.0 && e < 1.5) {
        if (strstr(desc->name, "SAT") || strstr(desc->name, "CNF"))
            return EQUIV_CLASS_CNFSAT;
    }
    return EQUIV_CLASS_COUNT;
}

int compare_equivalence_classes(equiv_class_id_t a, equiv_class_id_t b) {
    if (a == b) return 0;
    double da = class_diameter_exponent(a);
    double db = class_diameter_exponent(b);
    if (da >= 0 && db < 0) return 1;
    if (db >= 0 && da < 0) return -1;
    if (da < db) return 1;
    if (db < da) return -1;
    return 0;
}

void print_equivalence_class_summary(void) {
    ctx_init();
    printf("=== Fine-Grained Equivalence Classes ===\n");
    printf("Total registered problems: %d\n", ctx.num_problems);
    printf("Total registered reductions: %d\n", ctx.num_reductions);
    printf("\n");
    printf("%-3s %-30s %6s %6s %6s\n",
           "ID", "Class Name", "Thresh", "Best", "ConjLB");
    printf("--- ------------------------------ ------ ------ ------\n");
    for (int i = 0; i < EQUIV_CLASS_COUNT; i++) {
        if (!ctx.class_initialized[i]) continue;
        equiv_class_t *ec = &ctx.classes[i];
        printf("%2d  %-30s %6.2f %6.4f %6.2f\n",
               i, ec->canonical_name, ec->threshold_exponent,
               ec->best_algorithm_exponent, ec->conjectured_lower_bound);
    }
    printf("\n");
}

void print_class_reduction_graph(equiv_class_id_t class_id) {
    ctx_init();
    printf("=== Reduction Graph for %s ===\n", equiv_class_name(class_id));
    int32_t count = 0;
    for (int32_t i = 0; i < ctx.num_reductions; i++) {
        reduction_descriptor_t *r = &ctx.reductions[i];
        problem_id_t from = r->from_problem;
        if (ctx.problems[from].equiv_class != class_id) continue;
        printf("  %s --(%s, exp=%.2f, blow=%d)--> %s\n",
               ctx.problems[r->from_problem].name,
               r->reduction_type == FGR_SUBCUBIC ? "subcubic" :
               r->reduction_type == FGR_SUBQUADRATIC ? "subquad" : "3SUM",
               r->runtime_exponent, r->output_blowup,
               ctx.problems[r->to_problem].name);
        count++;
    }
    if (count == 0) printf("  (no reductions registered)\n");
    printf("  Total: %d reductions in this class\n", count);
}

/* ---- L8: Advanced Cross-Class Relations ---- */

bool class_a_at_least_as_hard_as_b(equiv_class_id_t a, equiv_class_id_t b,
                                    const char *assumption) {
    (void)assumption;
    /* If SETH holds, CNF-SAT >= OV-class >= 3SUM-class */
    if (a == EQUIV_CLASS_CNFSAT && b != EQUIV_CLASS_CNFSAT) return true;
    /* Subcubic and subquadratic are incomparable in general */
    return false;
}

double fg_landscape_distance(equiv_class_id_t a, equiv_class_id_t b) {
    problem_id_t ca = class_canonical_problem(a);
    problem_id_t cb = class_canonical_problem(b);
    if (ca == (problem_id_t)-1 || cb == (problem_id_t)-1) return -1.0;
    double dab = reduction_chain_exponent(ca, cb);
    double dba = reduction_chain_exponent(cb, ca);
    if (dab >= 0 && dba >= 0) return (dab < dba) ? dab : dba;
    if (dab >= 0) return dab;
    if (dba >= 0) return dba;
    return -1.0;
}

bool class_a_implies_class_b(equiv_class_id_t a, equiv_class_id_t b) {
    double dist = fg_landscape_distance(a, b);
    if (dist < 0) return false;
    double ta = ctx.classes[a].threshold_exponent;
    double tb = ctx.classes[b].threshold_exponent;
    return dist <= ta - tb + 0.5;
}