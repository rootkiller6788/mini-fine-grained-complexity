/* reduction_web.c -- Reduction Web Implementation */
#include "reduction_web.h"
#include "condlb.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

ExtendedReductionWeb* erw_create(int max_problems, int max_edges, int max_classes) {
    assert(max_problems > 0);
    ExtendedReductionWeb* erw = (ExtendedReductionWeb*)malloc(sizeof(ExtendedReductionWeb));
    if (!erw) return NULL;
    erw->problems = (ProblemNode*)calloc((size_t)max_problems, sizeof(ProblemNode));
    erw->n_problems = 0;
    erw->max_problems = max_problems;
    int me = max_edges > 0 ? max_edges : 16;
    erw->edges = (ReductionEdge*)malloc((size_t)me * sizeof(ReductionEdge));
    erw->n_edges = 0;
    erw->max_edges = me;
    erw->adjacency = (int**)calloc((size_t)max_problems, sizeof(int*));
    erw->adj_sizes = (int*)calloc((size_t)max_problems, sizeof(int));
    erw->adj_capacities = (int*)calloc((size_t)max_problems, sizeof(int));
    int mc = max_classes > 0 ? max_classes : 8;
    erw->classes = (EquivalenceClass*)malloc((size_t)mc * sizeof(EquivalenceClass));
    erw->n_classes = 0;
    erw->max_classes = mc;
    return erw;
}

void erw_free(ExtendedReductionWeb* erw) {
    if (!erw) return;
    for (int i = 0; i < erw->n_problems; i++) free(erw->adjacency[i]);
    for (int i = 0; i < erw->n_classes; i++) free(erw->classes[i].member_ids);
    free(erw->problems); free(erw->edges); free(erw->adjacency);
    free(erw->adj_sizes); free(erw->adj_capacities); free(erw->classes);
    free(erw);
}

static void erw_ensure_adj(ExtendedReductionWeb* erw, int from) {
    if (erw->adj_sizes[from] >= erw->adj_capacities[from]) {
        int nc = erw->adj_capacities[from] == 0 ? 4 : erw->adj_capacities[from] * 2;
        int* na = (int*)realloc(erw->adjacency[from], (size_t)nc * sizeof(int));
        if (!na) return;
        erw->adjacency[from] = na;
        erw->adj_capacities[from] = nc;
    }
}

int erw_add_problem(ExtendedReductionWeb* erw, const char* name, ProblemCategory cat,
                    double naive_exp, double best_known_exp) {
    assert(erw && name);
    if (erw->n_problems >= erw->max_problems) return -1;
    int id = erw->n_problems++;
    strncpy(erw->problems[id].name, name, 127);
    erw->problems[id].name[127] = 0;
    erw->problems[id].category = cat;
    erw->problems[id].naive_exponent = naive_exp;
    erw->problems[id].best_known_exponent = best_known_exp;
    erw->problems[id].conditional_lower_bound = 0.0;
    erw->problems[id].is_hard = 0;
    erw->problems[id].is_complete = 0;
    return id;
}

void erw_add_edge(ExtendedReductionWeb* erw, int from_id, int to_id,
                  double time_blowup, double polylog_blowup, int randomized, const char* paper_ref) {
    assert(erw);
    assert(from_id >= 0 && from_id < erw->n_problems);
    assert(to_id >= 0 && to_id < erw->n_problems);
    if (erw->n_edges >= erw->max_edges) {
        int nm = erw->max_edges * 2;
        ReductionEdge* ne = (ReductionEdge*)realloc(erw->edges, (size_t)nm * sizeof(ReductionEdge));
        if (!ne) return;
        erw->edges = ne;
        erw->max_edges = nm;
    }
    ReductionEdge* e = &erw->edges[erw->n_edges++];
    e->from_problem_id = from_id;
    e->to_problem_id = to_id;
    e->time_blowup = time_blowup;
    e->polylog_blowup = polylog_blowup;
    e->randomized = randomized;
    if (paper_ref) { strncpy(e->paper_ref, paper_ref, 255); e->paper_ref[255] = 0; }
    else e->paper_ref[0] = 0;
    erw_ensure_adj(erw, from_id);
    erw->adjacency[from_id][erw->adj_sizes[from_id]++] = to_id;
}

void erw_set_clb(ExtendedReductionWeb* erw, int problem_id, double lower_bound) {
    assert(erw);
    assert(problem_id >= 0 && problem_id < erw->n_problems);
    erw->problems[problem_id].conditional_lower_bound = lower_bound;
    erw->problems[problem_id].is_hard = (lower_bound > 0.0) ? 1 : 0;
}

void erw_mark_complete(ExtendedReductionWeb* erw, int problem_id) {
    assert(erw);
    assert(problem_id >= 0 && problem_id < erw->n_problems);
    erw->problems[problem_id].is_complete = 1;
}

/* Kosaraju SCC for equivalence classes */
static void dfs1(ExtendedReductionWeb* erw, int u, int* vis, int* ord, int* oi) {
    vis[u] = 1;
    for (int j = 0; j < erw->adj_sizes[u]; j++) {
        int v = erw->adjacency[u][j];
        if (!vis[v]) dfs1(erw, v, vis, ord, oi);
    }
    ord[(*oi)++] = u;
}

static void dfs2(ExtendedReductionWeb* erw, int u, int* vis, int* comp, int cid,
                 int** radj, int* rsz) {
    vis[u] = 1; comp[u] = cid;
    for (int j = 0; j < rsz[u]; j++) {
        int v = radj[u][j];
        if (!vis[v]) dfs2(erw, v, vis, comp, cid, radj, rsz);
    }
}

int erw_compute_classes(ExtendedReductionWeb* erw) {
    assert(erw);
    int n = erw->n_problems;
    if (n == 0) return 0;
    int* vis = (int*)calloc((size_t)n, sizeof(int));
    int* ord = (int*)malloc((size_t)n * sizeof(int));
    if (!vis || !ord) { free(vis); free(ord); return 0; }
    int oi = 0;
    for (int i = 0; i < n; i++) if (!vis[i]) dfs1(erw, i, vis, ord, &oi);
    int** radj = (int**)calloc((size_t)n, sizeof(int*));
    int* rsz = (int*)calloc((size_t)n, sizeof(int));
    int* rcaps = (int*)calloc((size_t)n, sizeof(int));
    if (!radj || !rsz || !rcaps) { free(vis); free(ord); free(radj); free(rsz); free(rcaps); return 0; }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < erw->adj_sizes[i]; j++) {
            int v = erw->adjacency[i][j];
            if (rsz[v] >= rcaps[v]) {
                int nc = rcaps[v] == 0 ? 4 : rcaps[v] * 2;
                int* nr = (int*)realloc(radj[v], (size_t)nc * sizeof(int));
                if (!nr) continue;
                radj[v] = nr; rcaps[v] = nc;
            }
            radj[v][rsz[v]++] = i;
        }
    }
    int* comp = (int*)malloc((size_t)n * sizeof(int));
    if (!comp) {
        for(int i=0;i<n;i++) free(radj[i]);
        free(radj); free(rsz); free(rcaps); free(vis); free(ord); return 0;
    }
    memset(vis, 0, (size_t)n * sizeof(int));
    int nc = 0;
    for (int i = n - 1; i >= 0; i--) {
        int u = ord[i];
        if (!vis[u]) { dfs2(erw, u, vis, comp, nc, radj, rsz); nc++; }
    }
    erw->n_classes = nc;
    if (erw->n_classes > erw->max_classes) {
        int nm = erw->n_classes;
        EquivalenceClass* ne = (EquivalenceClass*)realloc(erw->classes, (size_t)nm * sizeof(EquivalenceClass));
        if (ne) { erw->classes = ne; erw->max_classes = nm; }
    }
    for (int c = 0; c < nc; c++) {
        erw->classes[c].class_id = c;
        erw->classes[c].member_ids = NULL;
        erw->classes[c].n_members = 0;
        erw->classes[c].max_members = 0;
        erw->classes[c].hardness_exponent = 0.0;
        erw->classes[c].is_subcubic = 0;
        erw->classes[c].is_subquadratic = 0;
        erw->classes[c].representative[0] = 0;
    }
    for (int i = 0; i < n; i++) {
        int c = comp[i];
        EquivalenceClass* ec = &erw->classes[c];
        if (ec->n_members >= ec->max_members) {
            int nm = ec->max_members == 0 ? 4 : ec->max_members * 2;
            int* nmem = (int*)realloc(ec->member_ids, (size_t)nm * sizeof(int));
            if (!nmem) continue;
            ec->member_ids = nmem; ec->max_members = nm;
        }
        ec->member_ids[ec->n_members++] = i;
    }
    for (int c = 0; c < nc; c++) {
        if (erw->classes[c].n_members > 0) {
            int rep = erw->classes[c].member_ids[0];
            strncpy(erw->classes[c].representative, erw->problems[rep].name, 127);
            erw->classes[c].hardness_exponent = erw->problems[rep].naive_exponent;
        }
    }
    for (int i = 0; i < n; i++) free(radj[i]);
    free(radj); free(rsz); free(rcaps); free(comp); free(vis); free(ord);
    return nc;
}

int erw_find_class(const ExtendedReductionWeb* erw, int problem_id) {
    assert(erw);
    for (int c = 0; c < erw->n_classes; c++)
        for (int j = 0; j < erw->classes[c].n_members; j++)
            if (erw->classes[c].member_ids[j] == problem_id) return c;
    return -1;
}

int erw_class_members(const ExtendedReductionWeb* erw, int class_id, int* output, int max_output) {
    assert(erw && output);
    if (class_id < 0 || class_id >= erw->n_classes) return 0;
    int n = erw->classes[class_id].n_members;
    int count = n < max_output ? n : max_output;
    for (int i = 0; i < count; i++) output[i] = erw->classes[class_id].member_ids[i];
    return count;
}

/* Lower Bound Chain Construction */
LowerBoundChain* lbc_construct(const ExtendedReductionWeb* erw,
                               int hypothesis_problem_id,
                               int target_problem_id) {
    assert(erw);
    int n = erw->n_problems;
    int* parent = (int*)malloc((size_t)n * sizeof(int));
    int* edge_to = (int*)malloc((size_t)n * sizeof(int));
    int* visited = (int*)calloc((size_t)n, sizeof(int));
    int* queue = (int*)malloc((size_t)n * sizeof(int));
    if (!parent || !edge_to || !visited || !queue) {
        free(parent); free(edge_to); free(visited); free(queue); return NULL;
    }
    for (int i = 0; i < n; i++) { parent[i] = -1; edge_to[i] = -1; }
    int front = 0, back = 0;
    visited[hypothesis_problem_id] = 1;
    queue[back++] = hypothesis_problem_id;
    int found = 0;
    while (front < back && !found) {
        int u = queue[front++];
        if (u == target_problem_id) { found = 1; break; }
        for (int j = 0; j < erw->adj_sizes[u]; j++) {
            int v = erw->adjacency[u][j];
            if (!visited[v]) {
                visited[v] = 1;
                parent[v] = u;
                edge_to[v] = j;
                queue[back++] = v;
            }
        }
    }
    if (!found) { free(parent); free(edge_to); free(visited); free(queue); return NULL; }

    LowerBoundChain* chain = (LowerBoundChain*)malloc(sizeof(LowerBoundChain));
    if (!chain) { free(parent); free(edge_to); free(visited); free(queue); return NULL; }
    int* path = (int*)malloc((size_t)n * sizeof(int));
    int plen = 0;
    for (int v = target_problem_id; v != -1; v = parent[v]) path[plen++] = v;
    chain->length = plen - 1;
    chain->reductions = (WeightedReduction*)malloc((size_t)chain->length * sizeof(WeightedReduction));
    if (!chain->reductions) {
        free(chain); free(path); free(parent); free(edge_to); free(visited); free(queue); return NULL;
    }
    for (int k = 0; k < chain->length; k++) {
        int from = path[plen - 1 - k];
        int to   = path[plen - 2 - k];
        chain->reductions[k].from_id = from;
        chain->reductions[k].to_id = to;
        for (int e = 0; e < erw->n_edges; e++) {
            if (erw->edges[e].from_problem_id == from && erw->edges[e].to_problem_id == to) {
                chain->reductions[k].oracle_calls = 1.0;
                chain->reductions[k].input_blowup = erw->edges[e].time_blowup;
                chain->reductions[k].overhead_exponent = 0.0;
                break;
            }
        }
    }
    chain->final_exponent = lbc_compute_exponent(chain);
    chain->valid = lbc_validate(chain);
    free(path); free(parent); free(edge_to); free(visited); free(queue);
    return chain;
}

int lbc_validate(const LowerBoundChain* chain) {
    if (!chain) return 0;
    for (int i = 0; i < chain->length; i++) {
        if (chain->reductions[i].input_blowup <= 0.0) return 0;
    }
    return 1;
}

double lbc_compute_exponent(const LowerBoundChain* chain) {
    if (!chain || chain->length == 0) return 0.0;
    double exp = 2.0;
    for (int i = 0; i < chain->length; i++) {
        exp *= chain->reductions[i].input_blowup;
    }
    return exp;
}

void lbc_free(LowerBoundChain* chain) {
    if (!chain) return;
    free(chain->reductions);
    free(chain);
}

/* Subcubic and Subquadratic Equivalence */
int erw_same_subcubic_class(const ExtendedReductionWeb* erw, int p1, int p2) {
    assert(erw);
    int c1 = erw_find_class(erw, p1);
    int c2 = erw_find_class(erw, p2);
    if (c1 < 0 || c2 < 0) return 0;
    return (c1 == c2 && erw->classes[c1].is_subcubic) ? 1 : 0;
}

int erw_same_subquadratic_class(const ExtendedReductionWeb* erw, int p1, int p2) {
    assert(erw);
    int c1 = erw_find_class(erw, p1);
    int c2 = erw_find_class(erw, p2);
    if (c1 < 0 || c2 < 0) return 0;
    return (c1 == c2 && erw->classes[c1].is_subquadratic) ? 1 : 0;
}

void erw_propagate_hardness(ExtendedReductionWeb* erw) {
    assert(erw);
    for (int i = 0; i < erw->n_problems; i++) {
        if (erw->problems[i].is_hard) {
            for (int k = 0; k < erw->n_edges; k++) {
                if (erw->edges[k].from_problem_id == i) {
                    int to = erw->edges[k].to_problem_id;
                    double prop = erw->problems[i].conditional_lower_bound *
                                  erw->edges[k].time_blowup;
                    if (prop > erw->problems[to].conditional_lower_bound) {
                        erw->problems[to].conditional_lower_bound = prop;
                        erw->problems[to].is_hard = 1;
                    }
                }
            }
        }
    }
}

double erw_strongest_bound(const ExtendedReductionWeb* erw, int problem_id) {
    assert(erw);
    assert(problem_id >= 0 && problem_id < erw->n_problems);
    return erw->problems[problem_id].conditional_lower_bound;
}

/* Print functions */
void erw_print_subcubic_class(const ExtendedReductionWeb* erw) {
    (void)erw;
    printf("=== Subcubic Equivalence Class ===\n\n");
    printf("Problems requiring n^{3-o(1)} time under APSP conjecture:\n");
    printf("  - APSP\n");
    printf("  - Negative Triangle\n");
    printf("  - Graph Radius\n");
    printf("  - Median\n");
    printf("  - Betweenness Centrality\n");
    printf("  - Replacement Paths\n");
    printf("  - Second Shortest Path\n");
    printf("  - Minimum Weight Cycle\n");
    printf("Total known members: 32+ (Vassilevska-Williams-Williams, 2018)\n");
}

void erw_print_subquadratic_class(const ExtendedReductionWeb* erw) {
    (void)erw;
    printf("=== Subquadratic Equivalence Class ===\n\n");
    printf("Problems requiring n^{2-o(1)} time under SETH/OVC:\n");
    printf("  - Orthogonal Vectors\n");
    printf("  - Edit Distance\n");
    printf("  - LCS\n");
    printf("  - Frechet Distance\n");
    printf("  - Dynamic Time Warping\n");
    printf("  - Diameter\n");
    printf("  - Subgraph Isomorphism\n");
    printf("  - Hitting Set\n");
    printf("  - Regex Matching\n");
    printf("  - Longest Palindrome Subsequence\n");
    printf("Total known members: 20+ (Williams, 2018)\n");
}

void erw_print_summary(const ExtendedReductionWeb* erw) {
    assert(erw);
    printf("=== Extended Reduction Web Summary ===\n\n");
    printf("Problems: %d, Edges: %d, Classes: %d\n\n",
           erw->n_problems, erw->n_edges, erw->n_classes);
    printf("--- Problem Nodes ---\n");
    for (int i = 0; i < erw->n_problems; i++) {
        ProblemNode* p = &erw->problems[i];
        printf("[%d] %-25s naive=%.1f best=%.3f clb=%.1f hard=%d complete=%d\n",
               i, p->name, p->naive_exponent, p->best_known_exponent,
               p->conditional_lower_bound, p->is_hard, p->is_complete);
    }
    printf("\n--- Reduction Edges ---\n");
    for (int i = 0; i < erw->n_edges; i++) {
        ReductionEdge* e = &erw->edges[i];
        printf("[%d] %s -> %s (blowup=%.2f, plog=%.1f, rand=%d)\n",
               i, erw->problems[e->from_problem_id].name,
               erw->problems[e->to_problem_id].name,
               e->time_blowup, e->polylog_blowup, e->randomized);
        if (e->paper_ref[0])
            printf("    ref: %s\n", e->paper_ref);
    }
    printf("\n--- Equivalence Classes ---\n");
    for (int c = 0; c < erw->n_classes; c++) {
        EquivalenceClass* ec = &erw->classes[c];
        printf("Class %d: rep=%s members=%d exponent=%.1f\n",
               c, ec->representative, ec->n_members, ec->hardness_exponent);
        if (ec->n_members > 0) {
            printf("  Members: ");
            for (int j = 0; j < ec->n_members && j < 10; j++)
                printf("%s ", erw->problems[ec->member_ids[j]].name);
            if (ec->n_members > 10)
                printf("...");
            printf("\n");
        }
    }
}
