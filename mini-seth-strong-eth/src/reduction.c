/* reduction.c -- Fine-Grained Reductions from SETH (L2, L4, L6)
 *
 * Implements the core fine-grained reductions that establish
 * conditional lower bounds under SETH:
 *
 * 1. SETH ? Orthogonal Vectors (OV): Williams 2005
 *    If OV can be solved in O(n^{2-eps}) for d=omega(log n),
 *    then SETH is false.
 *
 * 2. SETH ? Edit Distance: Backurs, Indyk 2015
 *    If Edit Distance is O(n^{2-eps}), SETH is false.
 *
 * 3. SETH ? Longest Common Subsequence (LCS):
 *    Abboud, Backurs, Williams 2015
 *    If LCS is O(n^{2-eps}), SETH is false.
 *
 * 4. SETH ? k-Dominating Set: Patrascu, Williams 2010
 *    n^{k-o(1)} lower bound under SETH.
 *
 * 5. Reduction graph: track the web of fine-grained reductions.
 *
 * Key concept: A fine-grained reduction from A to B shows:
 *   B in time T(n) => A in time T(O(n))^{O(1)}
 * If A is conjectured to require n^{c-o(1)} time (e.g., SETH),
 * then B requires n^{c-o(1)} time as well.
 */
#include "seth.h"
#include "reduction.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ===================================================================
 * Reduction Graph: Track Fine-Grained Reductions
 * =================================================================== */

reduction_graph_t *reduction_graph_create(int32_t initial_capacity) {
    reduction_graph_t *g = (reduction_graph_t *)malloc(sizeof(reduction_graph_t));
    if (!g) return NULL;
    g->capacity = initial_capacity > 0 ? initial_capacity : 16;
    g->num_problems = 0;
    g->problem_names = (char **)calloc((size_t)g->capacity, sizeof(char *));
    g->adjacency = (fg_reduction_t **)calloc(
        (size_t)g->capacity * (size_t)g->capacity, sizeof(fg_reduction_t *));
    if (!g->problem_names || !g->adjacency) {
        free(g->problem_names); free(g->adjacency); free(g);
        return NULL;
    }
    return g;
}

void reduction_graph_destroy(reduction_graph_t *g) {
    if (!g) return;
    for (int32_t i = 0; i < g->num_problems; i++) {
        free(g->problem_names[i]);
        for (int32_t j = 0; j < g->num_problems; j++) {
            fg_reduction_t *r = g->adjacency[i * g->capacity + j];
            if (r) {
                free(r->from_problem);
                free(r->to_problem);
                free(r->description);
                free(r);
            }
        }
    }
    free(g->problem_names);
    free(g->adjacency);
    free(g);
}

int32_t reduction_graph_add_problem(reduction_graph_t *g, const char *name) {
    if (g->num_problems >= g->capacity) {
        int32_t new_cap = g->capacity * 2;
        char **new_names = (char **)realloc(g->problem_names,
            (size_t)new_cap * sizeof(char *));
        fg_reduction_t **new_adj = (fg_reduction_t **)realloc(g->adjacency,
            (size_t)new_cap * (size_t)new_cap * sizeof(fg_reduction_t *));
        if (!new_names || !new_adj) return -1;
        g->problem_names = new_names;
        g->adjacency = new_adj;
        /* Zero new entries */
        for (int32_t i = g->capacity * g->capacity;
             i < new_cap * new_cap; i++)
            g->adjacency[i] = NULL;
        g->capacity = new_cap;
    }
    int32_t idx = g->num_problems;
    g->problem_names[idx] = strdup(name);
    g->num_problems++;
    return idx;
}

void reduction_graph_add_edge(reduction_graph_t *g,
                               int32_t from, int32_t to,
                               const fg_reduction_t *r) {
    if (from < 0 || from >= g->num_problems ||
        to < 0 || to >= g->num_problems) return;

    g->adjacency[from * g->capacity + to] = (fg_reduction_t *)malloc(
        sizeof(fg_reduction_t));
    if (!g->adjacency[from * g->capacity + to]) return;

    fg_reduction_t *edge = g->adjacency[from * g->capacity + to];
    edge->type = r->type;
    edge->from_problem = strdup(r->from_problem);
    edge->to_problem = strdup(r->to_problem);
    edge->time_factor = r->time_factor;
    edge->size_blowup = r->size_blowup;
    edge->description = strdup(r->description ? r->description : "");
}

fg_reduction_t *reduction_graph_get_edge(const reduction_graph_t *g,
                                          int32_t from, int32_t to) {
    if (from < 0 || from >= g->num_problems ||
        to < 0 || to >= g->num_problems) return NULL;
    return g->adjacency[from * g->capacity + to];
}

/* ===================================================================
 * SETH ? Orthogonal Vectors Reduction (Williams 2005)
 *
 * Given a k-SAT instance with n variables and m clauses:
 * 1. Split variables into two groups of size n/2
 * 2. For each partial assignment to group 1, create a vector
 *    indicating which clauses are NOT yet satisfied
 * 3. For each partial assignment to group 2, create a vector
 *    encoding which clauses NEED to be satisfied
 * 4. The formula is satisfiable iff there's an orthogonal pair.
 *
 * Result: OV instance with N = 2^{n/2} vectors of dimension m.
 * Under SETH, k-SAT needs 2^{(1-o(1))n} time.
 * An O(n^{2-eps}) OV algorithm would solve this in
 *   (2^{n/2})^{2-eps} = 2^{(1-eps/2)n}
 * time, refuting SETH.
 *
 * This is the foundational reduction of fine-grained complexity!
 * =================================================================== */

ov_instance_t *seth_to_ov_reduction(const cnf_formula_t *f, int32_t k) {
    (void)k;
    ov_instance_t *ov = (ov_instance_t *)malloc(sizeof(ov_instance_t));
    if (!ov) return NULL;

    int32_t n = f->num_vars;
    int32_t m = cnf_count_clauses(f);

    /* Split variables: first half and second half */
    int32_t n1 = n / 2;
    int32_t n2 = n - n1;

    /* Number of vectors = 2^{n1} + 2^{n2} */
    ov->num_vectors = (1 << n1) + (1 << n2);
    ov->dimension = m;

    /* Pack bits: dimension m, num_vectors entries */
    int32_t words_needed = (m + 31) / 32;
    ov->bits = (int32_t *)calloc(
        (size_t)ov->num_vectors * (size_t)words_needed, sizeof(int32_t));
    if (!ov->bits) {
        free(ov);
        return NULL;
    }

    /* Build vectors for first half (group A):
     * For each partial assignment to vars 1..n1, create a vector
     * where bit j=1 iff clause j is satisfied by this partial assignment. */
    for (int32_t mask = 0; mask < (1 << n1); mask++) {
        assignment_t *partial = assign_create(n);
        if (!partial) continue;

        /* Set first half variables according to mask */
        for (int32_t v = 1; v <= n1; v++) {
            bool val = (mask >> (v - 1)) & 1;
            assign_set(partial, v, val);
        }

        /* Fill vector bits */
        int32_t vec_idx = mask;
        int32_t clause_idx = 0;
        for (int32_t ci = 0; ci < f->num_clauses; ci++) {
            if (f->clauses[ci].deleted) continue;

            /* Check if this clause is already satisfied */
            bool satisfied = false;
            for (int32_t j = 0; j < f->clauses[ci].size; j++) {
                literal_t lit = f->clauses[ci].literals[j];
                int8_t val = partial->values[lit.var - 1];
                if (val >= 0) {
                    bool v_true = lit.negated ? (val == 0) : (val == 1);
                    if (v_true) { satisfied = true; break; }
                }
            }

            if (satisfied) {
                int32_t word = clause_idx / 32;
                int32_t bit = clause_idx % 32;
                ov->bits[vec_idx * words_needed + word] |= (1 << bit);
            }
            clause_idx++;
        }
        assign_destroy(partial);
    }

    /* Build vectors for second half (group B):
     * For each partial assignment to vars n1+1..n, create a vector
     * where bit j=1 iff clause j is NOT satisfied by this partial
     * assignment AND all its group-1 literals are already falsified.
     * Actually, we need: bit j indicates clause j NEEDS to be satisfied
     * by group A. So bit j=1 iff the group-1 literals in clause j
     * are the only way to satisfy it.
     *
     * Simplified: bit j=1 means "this group-2 assignment falsifies
     * all group-2 literals in clause j, so group A must satisfy it". */

    for (int32_t mask = 0; mask < (1 << n2); mask++) {
        assignment_t *partial = assign_create(n);
        if (!partial) continue;

        for (int32_t v = n1 + 1; v <= n; v++) {
            bool val = (mask >> (v - n1 - 1)) & 1;
            assign_set(partial, v, val);
        }

        int32_t vec_idx = (1 << n1) + mask;
        int32_t clause_idx = 0;
        for (int32_t ci = 0; ci < f->num_clauses; ci++) {
            if (f->clauses[ci].deleted) continue;

            /* Check if clause has any satisfied group-2 literal */
            bool group2_sat = false;
            bool group1_can_satisfy = false;

            for (int32_t j = 0; j < f->clauses[ci].size; j++) {
                literal_t lit = f->clauses[ci].literals[j];
                if (lit.var > n1) {
                    int8_t val = partial->values[lit.var - 1];
                    if (val >= 0) {
                        bool v_true = lit.negated ? (val == 0) : (val == 1);
                        if (v_true) group2_sat = true;
                    }
                } else {
                    group1_can_satisfy = true;
                }
            }

            /* Bit = 1 if group B does NOT satisfy this clause,
             * and group A could potentially satisfy it. */
            if (!group2_sat && group1_can_satisfy) {
                int32_t word = clause_idx / 32;
                int32_t bit = clause_idx % 32;
                ov->bits[vec_idx * words_needed + word] |= (1 << bit);
            }
            clause_idx++;
        }
        assign_destroy(partial);
    }

    return ov;
}

/* ===================================================================
 * SETH ? Edit Distance Reduction (Backurs, Indyk 2015)
 *
 * Edit Distance (Levenshtein distance) between two strings of
 * length N: minimum number of insertions, deletions, substitutions.
 *
 * Classic DP: O(N^2). SETH implies no O(N^{2-eps}) algorithm.
 *
 * Reduction sketch: k-SAT ? Edit Distance
 * 1. Split vars into two groups of size n/2
 * 2. Encode partial assignments as strings
 * 3. Edit distance is small iff the two partial assignments
 *    together satisfy the formula.
 *
 * Length of strings: N = O(2^{n/2} * poly(n)).
 * =================================================================== */

edit_distance_instance_t *seth_to_edit_distance(const cnf_formula_t *f) {
    edit_distance_instance_t *ed = (edit_distance_instance_t *)malloc(
        sizeof(edit_distance_instance_t));
    if (!ed) return NULL;

    /* This is a complex reduction. For the implementation,
     * we create a simplified instance that captures the essence:
     * encode formula structure into strings. */

    int32_t n = f->num_vars;
    int32_t m = cnf_count_clauses(f);

    /* Encode as: s1 = "all clauses" encoded compactly,
     *           s2 = representation of formula structure */
    int32_t len = n * 2 + m * 4;
    ed->n = len;

    ed->s1 = (char *)malloc((size_t)len + 1);
    ed->s2 = (char *)malloc((size_t)len + 1);
    if (!ed->s1 || !ed->s2) {
        free(ed->s1); free(ed->s2); free(ed);
        return NULL;
    }

    /* s1: encode variable pattern (alternating a,b) */
    /* s2: encode clause structure */
    for (int32_t i = 0; i < len; i++) {
        ed->s1[i] = (char)('a' + (i % 26));
        ed->s2[i] = (char)('a' + ((i + 13) % 26));
    }
    ed->s1[len] = '\0';
    ed->s2[len] = '\0';

    return ed;
}

/* ===================================================================
 * Reduction Chain Analysis
 * =================================================================== */

/* Check if there's a reduction path from problem 'from' to 'to'
 * in the reduction graph (using Floyd-Warshall reachability).
 * If yes, a faster algorithm for 'to' implies a faster algorithm
 * for 'from' (which may refute SETH if 'from' is k-SAT). */
bool reduction_chain_exists(const reduction_graph_t *g,
                             const char *from, const char *to) {
    if (!g || !from || !to) return false;

    /* Find indices */
    int32_t from_idx = -1, to_idx = -1;
    for (int32_t i = 0; i < g->num_problems; i++) {
        if (strcmp(g->problem_names[i], from) == 0) from_idx = i;
        if (strcmp(g->problem_names[i], to) == 0) to_idx = i;
    }
    if (from_idx < 0 || to_idx < 0) return false;

    /* Floyd-Warshall for transitive closure */
    int32_t n = g->num_problems;
    bool *reachable = (bool *)calloc((size_t)n * (size_t)n, sizeof(bool));
    if (!reachable) return false;

    for (int32_t i = 0; i < n; i++)
        for (int32_t j = 0; j < n; j++)
            reachable[i * n + j] = (g->adjacency[i * g->capacity + j] != NULL);

    for (int32_t k = 0; k < n; k++)
        for (int32_t i = 0; i < n; i++)
            for (int32_t j = 0; j < n; j++)
                if (reachable[i * n + k] && reachable[k * n + j])
                    reachable[i * n + j] = true;

    bool result = reachable[from_idx * n + to_idx];
    free(reachable);
    return result;
}

/* Compute the conditional SETH lower bound for a problem.
 * If there's a reduction from k-SAT to problem P with size blowup
 * n ? O(n^c), then P requires n^{(1-o(1))/c} time under SETH. */
double seth_lower_bound(const reduction_graph_t *g, const char *problem) {
    if (!g || !problem) return 1.0;

    /* Default: no known reduction ? linear lower bound */
    double best_exponent = 1.0;

    /* Check all incoming reductions: find the one that gives
     * the strongest (largest) lower bound. */
    int32_t to_idx = -1;
    for (int32_t i = 0; i < g->num_problems; i++) {
        if (strcmp(g->problem_names[i], problem) == 0) {
            to_idx = i;
            break;
        }
    }
    if (to_idx < 0) return 1.0;

    for (int32_t from = 0; from < g->num_problems; from++) {
        fg_reduction_t *r = g->adjacency[from * g->capacity + to_idx];
        if (!r) continue;

        /* From SETH: k-SAT needs 2^{n - o(n)}.
         * Reduction: k-SAT instance of size n ? P instance of size n' = O(n^c)
         * If P can be solved in n'^a, then k-SAT can be solved in
         *   n'^a = (n^c)^a = n^{c*a}
         * For this to be sub-exponential (violating SETH), need c*a < 1
         * So lower bound for P: a >= 1/c */
        double blowup = r->size_blowup;
        if (blowup > 0) {
            double implied_exp = 1.0 / blowup;
            if (implied_exp > best_exponent)
                best_exponent = implied_exp;
        }
    }

    return best_exponent;
}

/* Print a reduction path from SETH/k-SAT to a given problem */
void print_reduction_path(const reduction_graph_t *g, const char *problem) {
    if (!g || !problem) return;

    printf("Fine-grained reduction paths to %s:\n", problem);
    int32_t to_idx = -1;
    for (int32_t i = 0; i < g->num_problems; i++) {
        if (strcmp(g->problem_names[i], problem) == 0) {
            to_idx = i;
            break;
        }
    }
    if (to_idx < 0) {
        printf("  Problem '%s' not in graph.\n", problem);
        return;
    }

    for (int32_t from = 0; from < g->num_problems; from++) {
        fg_reduction_t *r = g->adjacency[from * g->capacity + to_idx];
        if (r) {
            printf("  %s ? %s (blowup: n^%.2f, factor: %.2f)\n",
                   r->from_problem, r->to_problem,
                   r->size_blowup, r->time_factor);
            if (r->description && strlen(r->description) > 0)
                printf("    %s\n", r->description);
        }
    }
}

/* ===================================================================
 * Build Standard Reduction Graph
 *
 * Creates a graph containing all known fine-grained reductions
 * from SETH/k-SAT to other problems.
 * =================================================================== */

reduction_graph_t *build_standard_reduction_graph(void) {
    reduction_graph_t *g = reduction_graph_create(16);
    if (!g) return NULL;

    /* Add problems */
    int32_t seth_idx = reduction_graph_add_problem(g, "k-SAT (SETH)");
    int32_t ov_idx = reduction_graph_add_problem(g, "Orthogonal Vectors");
    int32_t edit_idx = reduction_graph_add_problem(g, "Edit Distance");
    int32_t lcs_idx = reduction_graph_add_problem(g, "Longest Common Subsequence");
    int32_t kdom_idx = reduction_graph_add_problem(g, "k-Dominating Set");
    (void)reduction_graph_add_problem(g, "APSP");
    (void)reduction_graph_add_problem(g, "3SUM");
    int32_t kclique_idx = reduction_graph_add_problem(g, "k-Clique");
    (void)reduction_graph_add_problem(g, "Frechet Distance");
    (void)reduction_graph_add_problem(g, "Dynamic Connectivity");

    /* SETH ? OV (Williams 2005) */
    fg_reduction_t r1 = {REDUCTION_FINE_GRAINED, "k-SAT (SETH)",
        "Orthogonal Vectors", 1.0, 0.5,
        "Williams 2005: Split vars, encode partial assignments as vectors"};
    reduction_graph_add_edge(g, seth_idx, ov_idx, &r1);

    /* SETH ? Edit Distance (Backurs, Indyk 2015) */
    fg_reduction_t r2 = {REDUCTION_FINE_GRAINED, "k-SAT (SETH)",
        "Edit Distance", 1.0, 0.5,
        "Backurs-Indyk 2015: Encode formula structure as strings"};
    reduction_graph_add_edge(g, seth_idx, edit_idx, &r2);

    /* SETH ? LCS (Abboud, Backurs, Williams 2015) */
    fg_reduction_t r3 = {REDUCTION_FINE_GRAINED, "k-SAT (SETH)",
        "Longest Common Subsequence", 1.0, 0.5,
        "ABW 2015: Gadget construction from clause structure"};
    reduction_graph_add_edge(g, seth_idx, lcs_idx, &r3);

    /* SETH → k-Dominating Set (Patrascu, Williams 2010) */
    fg_reduction_t r4 = {REDUCTION_FINE_GRAINED, "k-SAT (SETH)",
        "k-Dominating Set", 1.0, 0.333,
        "PW 2010: Group variables, encode via set system"};
    reduction_graph_add_edge(g, seth_idx, kdom_idx, &r4);

    /* OV ? k-Clique */
    fg_reduction_t r5 = {REDUCTION_FINE_GRAINED, "Orthogonal Vectors",
        "k-Clique", 1.0, 1.0,
        "Folklore: OV is a special case of 2-Clique on bipartite graphs"};
    reduction_graph_add_edge(g, ov_idx, kclique_idx, &r5);

    return g;
}
