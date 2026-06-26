/* condlb.c -- Conditional Lower Bounds: Core Implementation
 *
 * Implements the formal framework for conditional lower bounds (CLB).
 *
 * A conditional lower bound is a statement of the form:
 *   "If Hypothesis H is true, then Problem P requires time at least L(n)."
 *
 * This methodology, pioneered by Williams (2005) and systematized in
 * the 2010s, provides a way to prove meaningful lower bounds for
 * polynomial-time problems by conditioning on unproven (but widely
 * believed) hardness assumptions like SETH, the 3SUM conjecture,
 * and the APSP conjecture.
 *
 * The reduction web is the central data structure: a directed graph
 * where nodes are problems and edges are fine-grained reductions.
 * If problem A reduces to problem B, then the hardness of A "flows"
 * to B: a conditional lower bound for A becomes a conditional lower
 * bound for B.
 *
 * References:
 *   Williams (2018): "Fine-Grained Complexity" (textbook draft)
 *   Williams (2015): "Hardness of Easy Problems"
 *   Patrascu (2010): "Towards Polynomial Lower Bounds for Dynamic Problems"
 *   Abboud & Williams (2014): "Popular Conjectures Imply Strong Lower Bounds"
 *   Williams & Williams (2010): "Subcubic Equivalences"
 */

#include "condlb.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

/* ============================================================================
 * L1-L3: Hypothesis Database Management
 * ============================================================================ */

HypothesisDatabase* hyp_db_create(int max_hypotheses, int max_bounds) {
    assert(max_hypotheses > 0 && max_bounds > 0);
    HypothesisDatabase* db = (HypothesisDatabase*)malloc(sizeof(HypothesisDatabase));
    if (!db) return NULL;
    db->hypotheses = (Hypothesis*)calloc((size_t)max_hypotheses, sizeof(Hypothesis));
    if (!db->hypotheses) { free(db); return NULL; }
    db->n_hypotheses  = 0;
    db->max_hypotheses = max_hypotheses;
    db->bounds = (ConditionalLowerBound*)calloc((size_t)max_bounds, sizeof(ConditionalLowerBound));
    if (!db->bounds) { free(db->hypotheses); free(db); return NULL; }
    db->n_bounds   = 0;
    db->max_bounds = max_bounds;
    return db;
}

void hyp_db_free(HypothesisDatabase* db) {
    if (!db) return;
    if (db->hypotheses) {
        for (int i = 0; i < db->n_hypotheses; i++) {
            free(db->hypotheses[i].depends_on_ids);
            free(db->hypotheses[i].implies_ids);
        }
        free(db->hypotheses);
    }
    free(db->bounds);
    free(db);
}

int hyp_db_add_hypothesis(HypothesisDatabase* db,
                          const char* name,
                          const char* problem,
                          double exponent,
                          double eps_tol,
                          HypothesisStatus status,
                          int year,
                          const char* proposed_by) {
    assert(db != NULL && name != NULL && problem != NULL);
    if (db->n_hypotheses >= db->max_hypotheses) {
        int new_max = db->max_hypotheses * 2;
        Hypothesis* new_h = (Hypothesis*)realloc(db->hypotheses,
            (size_t)new_max * sizeof(Hypothesis));
        if (!new_h) return -1;
        memset(new_h + db->max_hypotheses, 0,
               (size_t)(new_max - db->max_hypotheses) * sizeof(Hypothesis));
        db->hypotheses = new_h;
        db->max_hypotheses = new_max;
    }
    int id = db->n_hypotheses++;
    Hypothesis* h = &db->hypotheses[id];
    strncpy(h->name, name, 127); h->name[127] = '\0';
    strncpy(h->problem, problem, 127); h->problem[127] = '\0';
    h->lower_bound_exponent = exponent;
    h->epsilon_tolerance    = eps_tol;
    h->status               = status;
    h->year_proposed         = year;
    strncpy(h->proposed_by, proposed_by, 127); h->proposed_by[127] = '\0';
    h->n_depends_on  = 0;
    h->depends_on_ids = NULL;
    h->n_implies     = 0;
    h->implies_ids   = NULL;
    return id;
}

void hyp_db_add_implication(HypothesisDatabase* db, int from_id, int to_id) {
    assert(db != NULL);
    assert(from_id >= 0 && from_id < db->n_hypotheses);
    assert(to_id >= 0 && to_id < db->n_hypotheses);

    /* Add to "implies" list of from_id */
    Hypothesis* from = &db->hypotheses[from_id];
    int* new_imp = (int*)realloc(from->implies_ids,
        (size_t)(from->n_implies + 1) * sizeof(int));
    if (!new_imp) return;
    new_imp[from->n_implies] = to_id;
    from->implies_ids = new_imp;
    from->n_implies++;

    /* Add to "depends_on" list of to_id */
    Hypothesis* to = &db->hypotheses[to_id];
    int* new_dep = (int*)realloc(to->depends_on_ids,
        (size_t)(to->n_depends_on + 1) * sizeof(int));
    if (!new_dep) return;
    new_dep[to->n_depends_on] = from_id;
    to->depends_on_ids = new_dep;
    to->n_depends_on++;
}

int hyp_db_add_bound(HypothesisDatabase* db,
                     int hyp_id,
                     const char* target_problem,
                     double exponent,
                     ReductionType rtype,
                     int tightness,
                     const char* reference) {
    assert(db != NULL && target_problem != NULL);
    if (db->n_bounds >= db->max_bounds) {
        int new_max = db->max_bounds * 2;
        ConditionalLowerBound* new_b = (ConditionalLowerBound*)realloc(
            db->bounds, (size_t)new_max * sizeof(ConditionalLowerBound));
        if (!new_b) return -1;
        memset(new_b + db->max_bounds, 0,
               (size_t)(new_max - db->max_bounds) * sizeof(ConditionalLowerBound));
        db->bounds = new_b;
        db->max_bounds = new_max;
    }
    int id = db->n_bounds++;
    ConditionalLowerBound* b = &db->bounds[id];
    b->hypothesis_id      = hyp_id;
    strncpy(b->target_problem, target_problem, 127);
    b->target_problem[127] = '\0';
    b->time_bound_exponent = exponent;
    b->reduction_type      = rtype;
    b->tightness           = tightness;
    b->proven              = 1;
    strncpy(b->reference, reference, 255);
    b->reference[255] = '\0';
    return id;
}

int hyp_db_find_hypothesis(const HypothesisDatabase* db, const char* name) {
    assert(db != NULL && name != NULL);
    for (int i = 0; i < db->n_hypotheses; i++) {
        if (strcmp(db->hypotheses[i].name, name) == 0) return i;
    }
    return -1;
}

void hyp_db_describe(const HypothesisDatabase* db, int hyp_id,
                     char* buf, int buf_sz) {
    assert(db != NULL && buf != NULL && buf_sz > 0);
    if (hyp_id < 0 || hyp_id >= db->n_hypotheses) {
        snprintf(buf, (size_t)buf_sz, "(invalid hypothesis id %d)", hyp_id);
        return;
    }
    const Hypothesis* h = &db->hypotheses[hyp_id];
    const char* status_str = "UNREFUTED";
    if (h->status == HYPOTHESIS_REFUTED)    status_str = "REFUTED";
    if (h->status == HYPOTHESIS_OPEN)       status_str = "OPEN";
    if (h->status == HYPOTHESIS_EQUIVALENT) status_str = "EQUIVALENT";
    snprintf(buf, (size_t)buf_sz,
             "%s (%s): %s req. n^{%.3f} time, eps-tol=%.6f [%d, %s, %s]",
             h->name, h->problem,
             h->lower_bound_exponent > 10 ?
             "EXPONENTIAL lower bound" : "POLYNOMIAL lower bound",
             h->lower_bound_exponent,
             h->epsilon_tolerance,
             h->year_proposed,
             h->proposed_by,
             status_str);
}

void hyp_db_print(const HypothesisDatabase* db) {
    assert(db != NULL);
    printf("=== Hypothesis Database ===\n");
    printf("Total hypotheses: %d, Total bounds: %d\n\n",
           db->n_hypotheses, db->n_bounds);

    printf("--- Hypotheses ---\n");
    for (int i = 0; i < db->n_hypotheses; i++) {
        const Hypothesis* h = &db->hypotheses[i];
        printf("[%d] %-12s | problem: %-14s | bound: n^{%.3f} | year: %d | by: %s\n",
               i, h->name, h->problem, h->lower_bound_exponent,
               h->year_proposed, h->proposed_by);
        if (h->n_depends_on > 0) {
            printf("     depends on: ");
            for (int j = 0; j < h->n_depends_on; j++) {
                printf("%s ", db->hypotheses[h->depends_on_ids[j]].name);
            }
            printf("\n");
        }
        if (h->n_implies > 0) {
            printf("     implies:    ");
            for (int j = 0; j < h->n_implies; j++) {
                printf("%s ", db->hypotheses[h->implies_ids[j]].name);
            }
            printf("\n");
        }
    }

    printf("\n--- Conditional Lower Bounds ---\n");
    for (int i = 0; i < db->n_bounds; i++) {
        const ConditionalLowerBound* b = &db->bounds[i];
        printf("[%d] Under %s: %s requires n^{%.3f} [%s, tight=%d]\n",
               i,
               db->hypotheses[b->hypothesis_id].name,
               b->target_problem,
               b->time_bound_exponent,
               b->reference,
               b->tightness);
    }
}

/* ============================================================================
 * L3: Reduction Web Implementation
 * ============================================================================ */

ReductionWeb* rw_create(int max_problems) {
    assert(max_problems > 0);
    ReductionWeb* rw = (ReductionWeb*)malloc(sizeof(ReductionWeb));
    if (!rw) return NULL;
    rw->problem_names = (char*)calloc((size_t)max_problems, 128);
    if (!rw->problem_names) { free(rw); return NULL; }
    rw->name_offsets = (int*)calloc((size_t)max_problems, sizeof(int));
    if (!rw->name_offsets) { free(rw->problem_names); free(rw); return NULL; }
    rw->n_problems   = 0;
    rw->max_problems = max_problems;
    rw->edges      = NULL;
    rw->n_edges    = 0;
    rw->max_edges  = 0;
    rw->adjacency    = (int**)calloc((size_t)max_problems, sizeof(int*));
    rw->adj_sizes    = (int*)calloc((size_t)max_problems, sizeof(int));
    rw->adj_capacities = (int*)calloc((size_t)max_problems, sizeof(int));
    if (!rw->adjacency || !rw->adj_sizes || !rw->adj_capacities) {
        free(rw->adjacency); free(rw->adj_sizes); free(rw->adj_capacities);
        free(rw->problem_names); free(rw->name_offsets); free(rw);
        return NULL;
    }
    return rw;
}

void rw_free(ReductionWeb* rw) {
    if (!rw) return;
    for (int i = 0; i < rw->n_problems; i++) {
        free(rw->adjacency[i]);
    }
    free(rw->adjacency);
    free(rw->adj_sizes);
    free(rw->adj_capacities);
    free(rw->edges);
    free(rw->problem_names);
    free(rw->name_offsets);
    free(rw);
}

static void rw_adj_add(ReductionWeb* rw, int from, int to) {
    if (rw->adj_sizes[from] >= rw->adj_capacities[from]) {
        int new_cap = (rw->adj_capacities[from] == 0) ? 4 : rw->adj_capacities[from] * 2;
        int* new_adj = (int*)realloc(rw->adjacency[from], (size_t)new_cap * sizeof(int));
        if (!new_adj) return;
        rw->adjacency[from] = new_adj;
        rw->adj_capacities[from] = new_cap;
    }
    rw->adjacency[from][rw->adj_sizes[from]++] = to;
}

int rw_add_problem(ReductionWeb* rw, const char* name) {
    assert(rw != NULL && name != NULL);
    if (rw->n_problems >= rw->max_problems) return -1;
    int id = rw->n_problems++;
    rw->name_offsets[id] = id * 128;
    strncpy(rw->problem_names + rw->name_offsets[id], name, 127);
    rw->problem_names[rw->name_offsets[id] + 127] = '\0';
    rw->adj_sizes[id] = 0;
    rw->adj_capacities[id] = 0;
    rw->adjacency[id] = NULL;
    return id;
}

void rw_add_reduction(ReductionWeb* rw,
                      int from_id, int to_id,
                      double time_blowup,
                      double polylog_blowup,
                      int randomized,
                      const char* paper_ref) {
    assert(rw != NULL);
    assert(from_id >= 0 && from_id < rw->n_problems);
    assert(to_id >= 0 && to_id < rw->n_problems);

    if (rw->n_edges >= rw->max_edges) {
        int new_max = (rw->max_edges == 0) ? 16 : rw->max_edges * 2;
        ReductionEdge* new_e = (ReductionEdge*)realloc(rw->edges,
            (size_t)new_max * sizeof(ReductionEdge));
        if (!new_e) return;
        rw->edges = new_e;
        rw->max_edges = new_max;
    }
    ReductionEdge* e = &rw->edges[rw->n_edges++];
    e->from_problem_id = from_id;
    e->to_problem_id   = to_id;
    e->time_blowup     = time_blowup;
    e->polylog_blowup  = polylog_blowup;
    e->randomized      = randomized;
    if (paper_ref) {
        strncpy(e->paper_ref, paper_ref, 255);
        e->paper_ref[255] = '\0';
    } else {
        e->paper_ref[0] = '\0';
    }
    rw_adj_add(rw, from_id, to_id);
}

/* Breadth-first search to find if a path exists from `from` to `to` */
int rw_path_exists(const ReductionWeb* rw, int from_id, int to_id) {
    assert(rw != NULL);
    if (from_id == to_id) return 1;
    int n = rw->n_problems;
    int* visited = (int*)calloc((size_t)n, sizeof(int));
    if (!visited) return 0;
    int* queue = (int*)malloc((size_t)n * sizeof(int));
    if (!queue) { free(visited); return 0; }
    int front = 0, back = 0;
    visited[from_id] = 1;
    queue[back++] = from_id;
    int found = 0;
    while (front < back) {
        int u = queue[front++];
        if (u == to_id) { found = 1; break; }
        for (int j = 0; j < rw->adj_sizes[u]; j++) {
            int v = rw->adjacency[u][j];
            if (!visited[v]) {
                visited[v] = 1;
                queue[back++] = v;
            }
        }
    }
    free(visited);
    free(queue);
    return found;
}

/* Compute transitive closure using Floyd-Warshall */
int** rw_transitive_closure(const ReductionWeb* rw, int* n_out) {
    assert(rw != NULL);
    int n = rw->n_problems;
    *n_out = n;
    int** closure = (int**)malloc((size_t)n * sizeof(int*));
    if (!closure) return NULL;
    for (int i = 0; i < n; i++) {
        closure[i] = (int*)calloc((size_t)n, sizeof(int));
        if (!closure[i]) {
            for (int k = 0; k < i; k++) free(closure[k]);
            free(closure);
            return NULL;
        }
        closure[i][i] = 1; /* reflexive */
        for (int j = 0; j < rw->adj_sizes[i]; j++) {
            closure[i][rw->adjacency[i][j]] = 1;
        }
    }
    /* Floyd-Warshall transitive closure */
    for (int k = 0; k < n; k++) {
        for (int i = 0; i < n; i++) {
            if (closure[i][k]) {
                for (int j = 0; j < n; j++) {
                    if (closure[k][j]) {
                        closure[i][j] = 1;
                    }
                }
            }
        }
    }
    return closure;
}

int rw_equivalence_class(const ReductionWeb* rw, int problem_id,
                         int* output, int max_output) {
    assert(rw != NULL && output != NULL);
    int n_problems = rw->n_problems;
    int** closure_in = (int**)malloc((size_t)n_problems * sizeof(int*));
    if (!closure_in) return 0;

    for (int i = 0; i < n_problems; i++) {
        closure_in[i] = (int*)calloc((size_t)n_problems, sizeof(int));
        if (!closure_in[i]) {
            for (int k = 0; k < i; k++) free(closure_in[k]);
            free(closure_in);
            return 0;
        }
    }

    /* Build reachability from each node via BFS */
    for (int start = 0; start < n_problems; start++) {
        int* visited = (int*)calloc((size_t)n_problems, sizeof(int));
        int* queue   = (int*)malloc((size_t)n_problems * sizeof(int));
        if (!visited || !queue) {
            free(visited); free(queue);
            for (int k = 0; k < n_problems; k++) free(closure_in[k]);
            free(closure_in);
            return 0;
        }
        int front = 0, back = 0;
        visited[start] = 1;
        queue[back++] = start;
        while (front < back) {
            int u = queue[front++];
            closure_in[start][u] = 1;
            for (int j = 0; j < rw->adj_sizes[u]; j++) {
                int v = rw->adjacency[u][j];
                if (!visited[v]) {
                    visited[v] = 1;
                    queue[back++] = v;
                }
            }
        }
        free(visited);
        free(queue);
    }

    int count = 0;
    for (int i = 0; i < n_problems && count < max_output; i++) {
        /* i is equivalent to problem_id iff they are mutually reachable */
        if (closure_in[problem_id][i] && closure_in[i][problem_id]) {
            output[count++] = i;
        }
    }

    for (int k = 0; k < n_problems; k++) free(closure_in[k]);
    free(closure_in);
    return count;
}

int rw_hardness_tier(const ReductionWeb* rw, int problem_id) {
    assert(rw != NULL);
    /* Find the longest path from any source (node with no incoming edges)
     * to this problem. The tier is the maximum depth. */
    int n = rw->n_problems;
    int* indeg = (int*)calloc((size_t)n, sizeof(int));
    if (!indeg) return -1;

    /* Compute in-degree of each node */
    for (int i = 0; i < rw->n_edges; i++) {
        indeg[rw->edges[i].to_problem_id]++;
    }

    /* Find all sources (in-degree 0) */
    int* depth = (int*)calloc((size_t)n, sizeof(int));
    int* queue = (int*)malloc((size_t)n * sizeof(int));
    if (!depth || !queue) {
        free(indeg); free(depth); free(queue); return -1;
    }
    int front = 0, back = 0;

    for (int i = 0; i < n; i++) {
        if (indeg[i] == 0) {
            depth[i] = 0;
            queue[back++] = i;
        } else {
            depth[i] = -1;
        }
    }

    while (front < back) {
        int u = queue[front++];
        for (int j = 0; j < rw->adj_sizes[u]; j++) {
            int v = rw->adjacency[u][j];
            if (depth[u] + 1 > depth[v]) {
                depth[v] = depth[u] + 1;
            }
            indeg[v]--;
            if (indeg[v] == 0) {
                queue[back++] = v;
            }
        }
    }

    int result = depth[problem_id];
    free(indeg); free(depth); free(queue);
    return result;
}

void rw_print(const ReductionWeb* rw) {
    assert(rw != NULL);
    printf("=== Reduction Web ===\n");
    printf("Problems: %d, Edges: %d\n\n", rw->n_problems, rw->n_edges);
    printf("--- Problem Nodes ---\n");
    for (int i = 0; i < rw->n_problems; i++) {
        const char* name = rw->problem_names + rw->name_offsets[i];
        printf("[%d] %s\n", i, name);
    }
    printf("\n--- Reduction Edges ---\n");
    for (int i = 0; i < rw->n_edges; i++) {
        const ReductionEdge* e = &rw->edges[i];
        const char* from_name = rw->problem_names + rw->name_offsets[e->from_problem_id];
        const char* to_name   = rw->problem_names + rw->name_offsets[e->to_problem_id];
        printf("[%d] %s -> %s (blowup=%.3f, plog=%.1f, rand=%d)\n",
               i, from_name, to_name,
               e->time_blowup, e->polylog_blowup, e->randomized);
        if (e->paper_ref[0]) {
            printf("    ref: %s\n", e->paper_ref);
        }
    }
}

/* ============================================================================
 * L3: Fine-Grained Time Bound Implementation
 * ============================================================================ */

FineGrainedTimeBound fgtb_make(double exponent, double polylog_power,
                                int is_upper, int is_exp, double exp_base) {
    FineGrainedTimeBound tb;
    tb.exponent         = exponent;
    tb.polylog_power    = polylog_power;
    tb.is_upper_bound   = is_upper;
    tb.is_exponential   = is_exp;
    tb.exponential_base = exp_base;
    return tb;
}

int fgtb_compare(FineGrainedTimeBound a, FineGrainedTimeBound b) {
    /* Exponential bounds: compare bases first */
    if (a.is_exponential && b.is_exponential) {
        if (a.exponential_base < b.exponential_base) return -1;
        if (a.exponential_base > b.exponential_base) return 1;
        /* same base, compare exponent */
        if (a.exponent < b.exponent) return -1;
        if (a.exponent > b.exponent) return 1;
        return 0;
    }
    if (a.is_exponential && !b.is_exponential) return 1;  /* exp > poly */
    if (!a.is_exponential && b.is_exponential) return -1;

    /* Both polynomial */
    if (a.exponent < b.exponent) return -1;
    if (a.exponent > b.exponent) return 1;
    if (a.polylog_power < b.polylog_power) return -1;
    if (a.polylog_power > b.polylog_power) return 1;
    return 0;
}

int fgtb_violates(FineGrainedTimeBound b, ConditionalLowerBound* lb) {
    if (!lb || !lb->proven) return 0;
    /* An upper bound b violates a lower bound lb if b's exponent
     * is strictly less than lb's exponent. */
    if (b.is_upper_bound && !b.is_exponential) {
        if (b.exponent < lb->time_bound_exponent) return 1;
    }
    return 0;
}

void fgtb_format(FineGrainedTimeBound tb, char* buf, int buf_sz) {
    if (tb.is_exponential) {
        if (tb.polylog_power > 0.0) {
            snprintf(buf, (size_t)buf_sz, "O(2^{%.4f n} * n^{%.2f})",
                     tb.exponential_base, tb.polylog_power);
        } else {
            snprintf(buf, (size_t)buf_sz, "O(2^{%.4f n})", tb.exponential_base);
        }
    } else {
        if (tb.polylog_power > 0.0) {
            snprintf(buf, (size_t)buf_sz, "O(n^{%.4f} * log^{%.1f} n)",
                     tb.exponent, tb.polylog_power);
        } else {
            snprintf(buf, (size_t)buf_sz, "O(n^{%.4f})", tb.exponent);
        }
    }
}

/* ============================================================================
 * L5: Standard Hypothesis Database Initialization
 * ============================================================================ */

/* Initialize the standard set of fine-grained hypotheses and their
 * relationships. This captures the known structure of the field. */
void hyp_db_init_standard(HypothesisDatabase* db) {
    assert(db != NULL);

    /* Core hypotheses */
    int eth_id = hyp_db_add_hypothesis(db, "ETH", "3-SAT",
        0.0, 0.0, HYPOTHESIS_UNREFUTED, 1999,
        "Impagliazzo & Paturi");
    /* ETH: 3-SAT requires time 2^{Omega(n)}. We encode this as
     * exponential_base > 0 meaning Omega(2^{delta n}). */

    int seth_id = hyp_db_add_hypothesis(db, "SETH", "k-SAT",
        0.0, 0.0, HYPOTHESIS_UNREFUTED, 2001,
        "Impagliazzo & Paturi");
    /* SETH: forall eps>0, exists k: k-SAT requires 2^{(1-eps)n} time. */

    int ovc_id = hyp_db_add_hypothesis(db, "OVC", "Orthogonal Vectors",
        2.0, 1e-9, HYPOTHESIS_UNREFUTED, 2005,
        "Williams");
    /* OVC: OV requires n^{2-o(1)} time. */

    int threesum_id = hyp_db_add_hypothesis(db, "3SUM", "3SUM",
        2.0, 1e-9, HYPOTHESIS_UNREFUTED, 2010,
        "Patrascu");
    /* 3SUM: requires n^{2-o(1)} time. */

    int apsp_id = hyp_db_add_hypothesis(db, "APSP", "All-Pairs Shortest Paths",
        3.0, 1e-9, HYPOTHESIS_UNREFUTED, 2010,
        "Williams & Williams");
    /* APSP: requires n^{3-o(1)} time. */

    int neth_id = hyp_db_add_hypothesis(db, "NSETH", "Nondeterministic k-SAT",
        0.0, 0.0, HYPOTHESIS_UNREFUTED, 2016,
        "Carmosino et al.");
    /* NSETH: nondeterministic SETH. */

    int eth_k_id = hyp_db_add_hypothesis(db, "ETH-k", "Exact k-SAT",
        1.0 / 3.0, 0.0, HYPOTHESIS_UNREFUTED, 2014,
        "Lokshtanov et al.");
    /* ETH for parameterized complexity: k-SAT with n vars, m clauses
     * requires time 2^{o(k)} * n^{O(1)}. */

    /* Implications: SETH => ETH */
    hyp_db_add_implication(db, seth_id, eth_id);

    /* SETH => OVC (Williams, 2005) */
    hyp_db_add_implication(db, seth_id, ovc_id);

    /* Add conditional lower bounds */
    hyp_db_add_bound(db, seth_id, "Edit Distance", 2.0,
        REDUCTION_FINE_GRAINED, 1,
        "Backurs & Indyk (STOC 2015)");

    hyp_db_add_bound(db, seth_id, "LCS", 2.0,
        REDUCTION_FINE_GRAINED, 1,
        "Abboud, Backurs & Williams (FOCS 2015)");

    hyp_db_add_bound(db, seth_id, "Frechet Distance", 2.0,
        REDUCTION_FINE_GRAINED, 1,
        "Bringmann (FOCS 2014)");

    hyp_db_add_bound(db, seth_id, "Dynamic Time Warping", 2.0,
        REDUCTION_FINE_GRAINED, 1,
        "Abboud, Backurs & Williams (FOCS 2015); Bringmann & Kunnemann (FOCS 2015)");

    hyp_db_add_bound(db, seth_id, "Diameter", 2.0,
        REDUCTION_FINE_GRAINED, 1,
        "Roditty & Williams (ICALP 2013)");

    hyp_db_add_bound(db, apsp_id, "Negative Triangle", 3.0,
        REDUCTION_FINE_GRAINED, 1,
        "Williams & Williams (FOCS 2010)");

    hyp_db_add_bound(db, apsp_id, "Graph Radius", 3.0,
        REDUCTION_FINE_GRAINED, 1,
        "Abboud, Grandoni, Williams (SODA 2015)");

    hyp_db_add_bound(db, apsp_id, "Betweenness Centrality", 3.0,
        REDUCTION_FINE_GRAINED, 1,
        "Abboud, Grandoni, Williams (SODA 2015)");

    hyp_db_add_bound(db, apsp_id, "Median", 3.0,
        REDUCTION_FINE_GRAINED, 1,
        "Abboud, Grandoni, Williams (SODA 2015)");

    hyp_db_add_bound(db, threesum_id, "3-Points-on-a-Line", 2.0,
        REDUCTION_FINE_GRAINED, 1,
        "Gajentaan & Overmars (1995)");

    hyp_db_add_bound(db, threesum_id, "Minimum Area Triangle", 2.0,
        REDUCTION_FINE_GRAINED, 1,
        "Gajentaan & Overmars (1995)");

    hyp_db_add_bound(db, ovc_id, "Bichromatic Closest Pair", 2.0,
        REDUCTION_FINE_GRAINED, 1,
        "Williams (2005)");
}

/* ============================================================================
 * L7: Applications - CLB Derivation for Dynamic Problems
 * ============================================================================ */

/* Compute the conditional lower bound for a dynamic graph problem.
 *
 * Many dynamic graph problems have conditional lower bounds from
 * the OMv (Online Matrix-Vector) conjecture (Henzinger et al., 2015)
 * or from SETH/OVC.
 *
 * This function, given a problem description, returns the strongest
 * known conditional lower bound exponent. */
double clb_dynamic_problem(const char* problem_name) {
    if (!problem_name) return 0.0;

    /* Known dynamic CLBs:
     *   - Dynamic Reachability: n^{2-o(1)} under SETH (Abboud-Williams, 2014)
     *   - Dynamic Bipartite Matching: n^{2-o(1)} under SETH
     *   - Dynamic Shortest Paths: n^{2-o(1)} under SETH
     *   - Dynamic Subgraph Connectivity: n^{2-o(1)} under OMv
     *   - Dynamic Maximum Flow: n^{2-o(1)} under SETH
     */

    if (strstr(problem_name, "Reachability") ||
        strstr(problem_name, "reachability")) return 2.0;
    if (strstr(problem_name, "Matching") ||
        strstr(problem_name, "matching")) return 2.0;
    if (strstr(problem_name, "Shortest Path") ||
        strstr(problem_name, "shortest path")) return 2.0;
    if (strstr(problem_name, "Subgraph") ||
        strstr(problem_name, "subgraph")) return 2.0;
    if (strstr(problem_name, "Max Flow") ||
        strstr(problem_name, "max flow")) return 2.0;
    if (strstr(problem_name, "Connectivity") ||
        strstr(problem_name, "connectivity")) return 2.0;

    return 0.0; /* unknown */
}

/* Check consistency of an algorithm with known CLBs.
 * Given an algorithm with exponent `algo_exp`, check whether it
 * would refute any known conditional lower bound.
 * Returns a bitmask of refuted hypotheses:
 *   bit 0: SETH
 *   bit 1: OVC
 *   bit 2: 3SUM
 *   bit 3: APSP
 */
int clb_check_consistency(const HypothesisDatabase* db,
                          const char* problem_name,
                          double algo_exp) {
    assert(db != NULL);
    int refuted = 0;

    for (int i = 0; i < db->n_bounds; i++) {
        const ConditionalLowerBound* b = &db->bounds[i];
        if (strcmp(b->target_problem, problem_name) == 0) {
            /* If algo_exp < bound, the algorithm would refute the hypothesis */
            if (algo_exp < b->time_bound_exponent - 1e-9) {
                int hyp_idx = b->hypothesis_id;
                const char* hname = db->hypotheses[hyp_idx].name;
                if (strcmp(hname, "SETH") == 0) refuted |= 1;
                if (strcmp(hname, "OVC") == 0)  refuted |= 2;
                if (strcmp(hname, "3SUM") == 0) refuted |= 4;
                if (strcmp(hname, "APSP") == 0) refuted |= 8;
            }
        }
    }
    return refuted;
}
