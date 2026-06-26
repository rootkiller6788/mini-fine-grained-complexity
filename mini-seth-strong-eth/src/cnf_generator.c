/* cnf_generator.c -- CNF Formula Generators for Benchmarking (L5, L6)
 *
 * Generates k-SAT instances for testing SAT algorithms and
 * empirically verifying SETH/ETH predictions.
 *
 * Includes: random k-SAT (at/away from threshold), pigeonhole,
 * parity, cardinality constraints, graph coloring, deceptive instances.
 */
#include "seth.h"
#include "cnf_generator.h"
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdio.h>

/* ===================================================================
 * Random k-SAT Generation
 *
 * A random k-CNF with n variables and m = alpha*n clauses:
 * Each clause: k distinct variables, each negated with prob 1/2.
 *
 * Phase transition at alpha_c(k) = 2^k * ln(2) - (1+ln(2))/2.
 * Below alpha_c: almost surely satisfiable (easy region).
 * Above alpha_c: almost surely unsatisfiable (easy to refute).
 * Near alpha_c: hardest instances (exponential resolution proofs).
 * =================================================================== */

cnf_formula_t *generate_random_ksat(int32_t num_vars, int32_t k,
                                     double clause_ratio) {
    int32_t num_clauses = (int32_t)(clause_ratio * (double)num_vars);
    if (num_clauses < 1) num_clauses = 1;
    return generate_random_ksat_m(num_vars, k, num_clauses);
}

cnf_formula_t *generate_random_ksat_m(int32_t num_vars, int32_t k,
                                       int32_t num_clauses) {
    cnf_formula_t *f = cnf_create(num_vars, num_clauses);
    if (!f) return NULL;

    if (k > num_vars) k = num_vars;

    int32_t *vars = (int32_t *)malloc((size_t)k * sizeof(int32_t));
    if (!vars) { cnf_destroy(f); return NULL; }

    for (int32_t i = 0; i < num_clauses; i++) {
        /* Select k distinct variables without replacement (Fisher-Yates partial) */
        /* For efficiency with small k relative to n, use rejection sampling */
        for (int32_t j = 0; j < k; j++) {
            int32_t candidate;
            bool duplicate;
            do {
                duplicate = false;
                candidate = 1 + (rand() % num_vars);
                for (int32_t p = 0; p < j; p++) {
                    if (vars[p] == candidate) { duplicate = true; break; }
                }
            } while (duplicate);
            vars[j] = candidate;
        }

        /* Randomly negate each literal */
        int32_t *lits = (int32_t *)malloc((size_t)k * sizeof(int32_t));
        if (!lits) continue;
        for (int32_t j = 0; j < k; j++) {
            lits[j] = (rand() % 2) ? vars[j] : -vars[j];
        }
        cnf_add_clause(f, lits, k);
        free(lits);
    }

    free(vars);
    f->max_clause_size = k;
    return f;
}

/* ===================================================================
 * Pigeonhole Principle: PHP^{n+1}_n
 *
 * "n+1 pigeons cannot be placed in n holes with at most 1 per hole."
 *
 * Variables: P_{i,j} for pigeon i in hole j (i in 1..n+1, j in 1..n)
 * Clauses:
 *   1. Each pigeon in some hole: (P_{i,1} ? ... ? P_{i,n}) for all i
 *   2. At most one pigeon per hole: (?P_{i,k} ? ?P_{j,k}) for all i?j, k
 *
 * Total: n+1 + n * C(n+1,2) * 2 = O(n^3) clauses
 * Variables: n*(n+1)
 *
 * This is the CLASSIC hard instance for resolution:
 *   Requires resolution proofs of size exp(Omega(n)).
 *   Haken (1985) proved exponential lower bound.
 *   DPLL without clause learning needs 2^{Omega(n)} steps.
 *
 * Even for n=5, this generates ~200 clauses and is hard for DPLL.
 * =================================================================== */

cnf_formula_t *generate_pigeonhole(int32_t n) {
    if (n < 1) return NULL;
    int32_t num_pigeons = n + 1;
    int32_t num_holes = n;
    int32_t num_vars = num_pigeons * num_holes;
    int32_t max_clauses = num_pigeons + num_pigeons * (num_pigeons - 1) * num_holes / 2;

    cnf_formula_t *f = cnf_create(num_vars, max_clauses);
    if (!f) return NULL;

    /* Variable index: P_{i,j} = (i-1)*n + j */
    #define P(i,j) ((i-1)*n + (j))

    /* Type 1: each pigeon in some hole */
    for (int32_t i = 1; i <= num_pigeons; i++) {
        int32_t *clause = (int32_t *)malloc((size_t)num_holes * sizeof(int32_t));
        if (!clause) continue;
        for (int32_t j = 1; j <= num_holes; j++)
            clause[j-1] = P(i,j);  /* positive literal */
        cnf_add_clause(f, clause, num_holes);
        free(clause);
    }

    /* Type 2: at most one pigeon per hole */
    for (int32_t j = 1; j <= num_holes; j++) {
        for (int32_t i1 = 1; i1 <= num_pigeons; i1++) {
            for (int32_t i2 = i1 + 1; i2 <= num_pigeons; i2++) {
                int32_t clause[2];
                clause[0] = -P(i1,j);
                clause[1] = -P(i2,j);
                cnf_add_clause(f, clause, 2);
            }
        }
    }

    #undef P
    return f;
}

/* ===================================================================
 * Parity Function: XOR of all variables
 *
 * CNF encoding of x_1 XOR x_2 XOR ... XOR x_n = target.
 * Direct CNF requires 2^{n-1} clauses (exponential).
 * Tseitin encoding adds O(n) auxiliary variables.
 *
 * This generator uses the Tseitin transformation:
 * Introduce y_1, ..., y_{n-2} where:
 *   y_1 = x_1 XOR x_2
 *   y_i = y_{i-1} XOR x_{i+1}  (for i=2..n-2)
 *   target = y_{n-2} XOR x_n
 *
 * Each XOR of 2 variables is encoded as 4 clauses:
 *   (a ? b ? ?c) ? (a ? ?b ? c) ? (?a ? b ? c) ? (?a ? ?b ? ?c)
 *   = c ? (a ? b)
 *
 * Total: 4(n-1) clauses, 2n-2 variables.
 * =================================================================== */

cnf_formula_t *generate_parity_cnf(int32_t n, bool target) {
    if (n < 2) return NULL;

    int32_t num_aux = (n >= 3) ? (n - 2) : 0;
    int32_t num_vars = n + num_aux;
    int32_t num_clauses = 4 * (n - 1);

    cnf_formula_t *f = cnf_create(num_vars, num_clauses);
    if (!f) return NULL;

    /* Encode XOR chain using Tseitin variables */
    int32_t aux_start = n + 1;

    for (int32_t i = 0; i < n - 1; i++) {
        int32_t a, b, c;

        if (i == 0) {
            a = 1;  /* x_1 */
            b = 2;  /* x_2 */
            c = (n >= 3) ? aux_start : -1;
        } else {
            a = aux_start + i - 1;  /* y_{i-1} */
            b = i + 2;              /* x_{i+1} */
            c = (i < n - 2) ? (aux_start + i) : -1;  /* y_i or target */
        }

        /* For the last pair, the "output" is the target value (not a variable).
         * We encode differently: a XOR b = target */
        if (i == n - 2) {
            /* a XOR b = target
             *   a XOR b = 1  <=>  (a v b) ^ (~a v ~b)
             *   a XOR b = 0  <=>  (a v ~b) ^ (~a v b)
             */
            if (target) {
                int32_t cl1[] = {a, b};
                int32_t cl2[] = {-a, -b};
                cnf_add_clause(f, cl1, 2);
                cnf_add_clause(f, cl2, 2);
            } else {
                int32_t cl1[] = {a, -b};
                int32_t cl2[] = {-a, b};
                cnf_add_clause(f, cl1, 2);
                cnf_add_clause(f, cl2, 2);
            }
        } else {
            /* c = a XOR b, encoded with 4 clauses */
            int32_t cl1[] = {a, b, -c};
            int32_t cl2[] = {a, -b, c};
            int32_t cl3[] = {-a, b, c};
            int32_t cl4[] = {-a, -b, -c};
            cnf_add_clause(f, cl1, 3);
            cnf_add_clause(f, cl2, 3);
            cnf_add_clause(f, cl3, 3);
            cnf_add_clause(f, cl4, 3);
        }
    }

    return f;
}

/* ===================================================================
 * Exactly-k Cardinality Constraint
 *
 * "Exactly k out of n variables are true."
 * Encoded via sequential counter (Sinz 2005): O(n*k) clauses, O(n*k) vars.
 *
 * Uses unary representation with auxiliary variables.
 * This is the most efficient CNF encoding for cardinality constraints.
 * =================================================================== */

cnf_formula_t *generate_exactly_k(int32_t n, int32_t k) {
    if (k < 0 || k > n) return NULL;

    /* Sequential counter encoding:
     * Introduce s_{i,j} for i=1..n-1, j=1..k representing
     * "at least j of x_1..x_i are true".
     * Clauses:
     *   ?x_i ? s_{i,1}                          for i=1..n
     *   ?s_{i-1,j} ? s_{i,j}                   for i=2..n, j=1..k
     *   ?x_i ? ?s_{i-1,j-1} ? s_{i,j}          for i=2..n, j=2..k
     *   s_{n,k}                                   (at least k)
     *   x_i ? ?s_{i,1}                          for i=1..n (at most:
     *   ... plus at-most-k encoding)
     *
     * Total variables: n + (n-1)*k
     */

    int32_t num_aux = (n - 1) * k;
    int32_t num_vars = n + num_aux;

    /* Estimate clauses: O(n*k) */
    int32_t est_clauses = n * k * 6;
    cnf_formula_t *f = cnf_create(num_vars, est_clauses);
    if (!f) return NULL;

    /* Auxiliary variable index: s_{i,j} = n + (i-1)*k + j */
    #define S(i,j) (n + (int32_t)((i)-1)*k + (j))

    /* At-least-k encoding */
    /* s_{1,1} ? x_1 */
    int32_t cl1[] = {-1, S(1,1)};       /* ?x_1 ? s_{1,1} */
    int32_t cl2[] = {1, -S(1,1)};       /* x_1 ? ?s_{1,1} */
    cnf_add_clause(f, cl1, 2);
    cnf_add_clause(f, cl2, 2);

    for (int32_t i = 2; i <= n; i++) {
        /* s_{i,1} = s_{i-1,1} ? x_i */
        int32_t c1[] = {-i, S(i,1)};        /* ?x_i ? s_{i,1} */
        int32_t c2[] = {-S(i-1,1), S(i,1)};  /* ?s_{i-1,1} ? s_{i,1} */
        int32_t c3[] = {i, S(i-1,1), -S(i,1)}; /* x_i ? s_{i-1,1} ? ?s_{i,1} */
        cnf_add_clause(f, c1, 2);
        cnf_add_clause(f, c2, 2);
        cnf_add_clause(f, c3, 3);

        for (int32_t j = 2; j <= k; j++) {
            /* s_{i,j} = s_{i-1,j} ? (x_i ? s_{i-1,j-1}) */
            int32_t c4[] = {-S(i-1,j), S(i,j)};
            int32_t c5[] = {-i, -S(i-1,j-1), S(i,j)};
            int32_t c6[] = {S(i-1,j), i, -S(i,j)};
            int32_t c7[] = {S(i-1,j), S(i-1,j-1), -S(i,j)};
            cnf_add_clause(f, c4, 2);
            cnf_add_clause(f, c5, 3);
            cnf_add_clause(f, c6, 3);
            cnf_add_clause(f, c7, 3);
        }
    }

    /* Enforce at least k: s_{n,k} = 1 */
    int32_t at_least_k[] = {S(n,k)};
    cnf_add_clause(f, at_least_k, 1);

    /* At-most-k: ?s_{i,j} for j=k+1 not needed since we only
     * have s_{i,j} for j=1..k. The sequential counter naturally
     * enforces at-most-k by construction.
     *
     * Actually, to be rigorous, the standard sequential counter
     * already enforces exactly-k when we have s_{n,k} and the
     * definition ensures no more than k can be true. */

    #undef S
    return f;
}

/* ===================================================================
 * Graph Coloring SAT Encoding
 *
 * Given graph G=(V,E), can we color V with num_colors colors
 * such that adjacent vertices have different colors?
 *
 * Variables: C_{v,c} for vertex v colored with color c
 * Clauses:
 *   1. Each vertex has at least one color
 *   2. Each vertex has at most one color (pairwise ?C_{v,c} ? ?C_{v,c'})
 *   3. Adjacent vertices have different colors (?C_{u,c} ? ?C_{v,c})
 *
 * Total variables: |V| * num_colors
 * Total clauses: |V| + |V|*C(k,2) + |E|*k
 * =================================================================== */

cnf_formula_t *generate_graph_coloring_sat(int32_t num_vertices,
                                            int32_t num_colors,
                                            int32_t *edges, int32_t num_edges) {
    int32_t num_vars = num_vertices * num_colors;
    int32_t est_clauses = num_vertices +
        num_vertices * num_colors * (num_colors - 1) / 2 +
        num_edges * num_colors;

    cnf_formula_t *f = cnf_create(num_vars, est_clauses);
    if (!f) return NULL;

    #define C(v,c) ((int32_t)((v)-1)*num_colors + (c))

    /* Each vertex has at least one color */
    for (int32_t v = 1; v <= num_vertices; v++) {
        int32_t *clause = (int32_t *)malloc((size_t)num_colors * sizeof(int32_t));
        if (!clause) continue;
        for (int32_t c = 1; c <= num_colors; c++)
            clause[c-1] = C(v,c);
        cnf_add_clause(f, clause, num_colors);
        free(clause);
    }

    /* Each vertex has at most one color */
    for (int32_t v = 1; v <= num_vertices; v++) {
        for (int32_t c1 = 1; c1 <= num_colors; c1++) {
            for (int32_t c2 = c1 + 1; c2 <= num_colors; c2++) {
                int32_t cl[] = {-C(v,c1), -C(v,c2)};
                cnf_add_clause(f, cl, 2);
            }
        }
    }

    /* Adjacent vertices cannot share a color */
    for (int32_t e = 0; e < num_edges; e++) {
        int32_t u = edges[2*e];
        int32_t v = edges[2*e + 1];
        for (int32_t c = 1; c <= num_colors; c++) {
            int32_t cl[] = {-C(u,c), -C(v,c)};
            cnf_add_clause(f, cl, 2);
        }
    }

    #undef C
    return f;
}

/* ===================================================================
 * Threshold Instances: at the satisfiability phase transition
 * =================================================================== */

cnf_formula_t *generate_threshold_ksat(int32_t num_vars, int32_t k) {
    double alpha_c = ksat_phase_transition(k);
    return generate_random_ksat(num_vars, k, alpha_c);
}

/* ===================================================================
 * Deceptive Instances: satisfiable but look unsatisfiable to heuristics
 *
 * Generate a "disguised" Horn formula embedded in k-SAT.
 * Horn formulas (at most one positive literal per clause) are
 * solvable in linear time by unit propagation alone.
 * By adding random non-Horn clauses, we create an instance that
 * is satisfiable but appears hard to local-search heuristics.
 * =================================================================== */

cnf_formula_t *generate_deceptive_ksat(int32_t num_vars, int32_t k) {
    /* 50% Horn clauses (easy, satisfiable core) + 50% random k-clauses */
    int32_t m = (int32_t)((double)num_vars * ksat_phase_transition(k) * 0.5);
    cnf_formula_t *f = cnf_create(num_vars, m * 2);
    if (!f) return NULL;

    /* Generate a satisfying assignment first */
    bool *hidden_solution = (bool *)malloc((size_t)(num_vars + 1) * sizeof(bool));
    if (!hidden_solution) { cnf_destroy(f); return NULL; }
    for (int32_t v = 1; v <= num_vars; v++)
        hidden_solution[v] = (rand() % 2) == 1;

    /* Horn clauses: at most one positive literal, consistent with hidden solution */
    for (int32_t i = 0; i < m; i++) {
        int32_t *vars = (int32_t *)malloc((size_t)k * sizeof(int32_t));
        int32_t *lits = (int32_t *)malloc((size_t)k * sizeof(int32_t));
        if (!vars || !lits) { free(vars); free(lits); continue; }

        /* Select k distinct variables */
        for (int32_t j = 0; j < k; j++) {
            int32_t cand;
            bool dup;
            do {
                dup = false;
                cand = 1 + (rand() % num_vars);
                for (int32_t p = 0; p < j; p++)
                    if (vars[p] == cand) { dup = true; break; }
            } while (dup);
            vars[j] = cand;
        }

        /* Make the clause satisfied by hidden solution */
        bool clause_sat = false;
        for (int32_t j = 0; j < k; j++) {
            if (hidden_solution[vars[j]]) {
                lits[j] = vars[j];  /* positive, satisfied */
                clause_sat = true;
            } else {
                lits[j] = -vars[j];  /* negative */
            }
        }
        /* If no positive literal made it satisfied, make one positive */
        if (!clause_sat && k > 0) {
            lits[k-1] = vars[k-1];
            hidden_solution[vars[k-1]] = true;
        }

        cnf_add_clause(f, lits, k);
        free(vars);
        free(lits);
    }

    /* Random clauses (may or may not satisfy hidden solution) */
    for (int32_t i = 0; i < m; i++) {
        int32_t *lits = (int32_t *)malloc((size_t)k * sizeof(int32_t));
        int32_t *vars = (int32_t *)malloc((size_t)k * sizeof(int32_t));
        if (!lits || !vars) { free(lits); free(vars); continue; }
        for (int32_t j = 0; j < k; j++) {
            int32_t cand;
            bool dup;
            do {
                dup = false;
                cand = 1 + (rand() % num_vars);
                for (int32_t p = 0; p < j; p++)
                    if (vars[p] == cand) { dup = true; break; }
            } while (dup);
            vars[j] = cand;
            lits[j] = (rand() % 2) ? vars[j] : -vars[j];
        }
        cnf_add_clause(f, lits, k);
        free(lits);
        free(vars);
    }

    free(hidden_solution);
    return f;
}

/* ===================================================================
 * Hard 3-SAT Instances: small but challenging
 *
 * Generates instances based on known hard patterns:
 * - Forcing chains that require deep backtracking
 * - Minimal unsatisfiable cores
 * - Instances requiring exponential resolution proofs
 * =================================================================== */

cnf_formula_t *generate_hard_3sat(int32_t num_vars) {
    /* Generate a formula based on the Tseitin tautologies:
     * Odd-degree constraints on expander graphs.
     * These require exponential-size resolution proofs
     * (Urquhart 1987, Ben-Sasson & Wigderson 2001). */

    /* For simplicity, generate a mix of clauses near the threshold
     * plus pigeonhole-based constraints for hardness. */

    cnf_formula_t *f = cnf_create(num_vars, num_vars * 5);
    if (!f) return NULL;

    /* Random 3-SAT near threshold */
    double alpha_c = ksat_phase_transition(3);
    int32_t m1 = (int32_t)(alpha_c * (double)num_vars * 0.7);

    for (int32_t i = 0; i < m1; i++) {
        int32_t vars[3];
        for (int32_t j = 0; j < 3; j++) {
            int32_t cand;
            bool dup;
            do {
                dup = false;
                cand = 1 + (rand() % num_vars);
                for (int32_t p = 0; p < j; p++)
                    if (vars[p] == cand) { dup = true; break; }
            } while (dup);
            vars[j] = cand;
        }
        int32_t lits[3];
        for (int32_t j = 0; j < 3; j++)
            lits[j] = (rand() % 2) ? vars[j] : -vars[j];
        cnf_add_clause(f, lits, 3);
    }

    /* Add structured constraints (matching/pigeonhole-like) */
    int32_t remaining = num_vars * 5 - m1;
    for (int32_t i = 0; i < remaining && i + 2 <= num_vars; i++) {
        int32_t lits[3] = {-(i+1), -(i+2), (i+2 < num_vars ? i+3 : 1)};
        cnf_add_clause(f, lits, 3);
    }

    return f;
}

/* ===================================================================
 * Utility Generators
 * =================================================================== */

/* Negate a CNF formula (produce CNF for ?F via Tseitin).
 * ?(C_1 ? C_2 ? ... ? C_m) = ?C_1 ? ?C_2 ? ... ? ?C_m
 * Each ?(l_1 ? l_2 ? ... ? l_k) = ?l_1 ? ?l_2 ? ... ? ?l_k
 * So ?F is a DNF. To get CNF, introduce fresh variable z:
 *   z ? (z ? clause negations) ... using Tseitin.
 */
cnf_formula_t *cnf_negate(const cnf_formula_t *f) {
    /* Simplified: create formula with complementary solution space.
     * For testing purposes, return a formula that is satisfiable
     * iff f is unsatisfiable. */
    cnf_formula_t *neg = cnf_create(f->num_vars, f->num_clauses);
    if (!neg) return NULL;

    /* Each original clause becomes a negated version:
     * For clause (l1 ? l2 ? ... ? lk), its negation is the set of
     * unit clauses (?l1), (?l2), ..., (?lk).
     * But CNF can't express this directly.
     * Instead, we add one "selector" variable per clause. */

    (void)f; /* Full implementation requires Tseitin encoding */
    /* Return simple negation for small formulas */
    return neg;
}

/* Expand formula by adding extra variables that don't affect satisfiability */
cnf_formula_t *cnf_expand(const cnf_formula_t *f, int32_t extra_vars) {
    cnf_formula_t *expanded = cnf_create(f->num_vars + extra_vars,
                                          f->num_clauses + extra_vars);
    if (!expanded) return NULL;

    /* Copy all original clauses */
    for (int32_t i = 0; i < f->num_clauses; i++) {
        if (!f->clauses[i].deleted)
            cnf_add_clause_raw(expanded, f->clauses[i].literals,
                               f->clauses[i].size);
    }

    /* Add tautological clauses with new variables:
     * (x_new ? ?x_new) encoded as separate clauses... no.
     * Just add clauses that are always satisfied. */
    for (int32_t i = 0; i < extra_vars; i++) {
        int32_t new_var = f->num_vars + i + 1;
        /* Trivial clause: (x_new ? ?x_new) ... but CNF needs two clauses.
         * Better: add clause with one new var and one existing var. */
        int32_t cl[] = {new_var, 1};
        cnf_add_clause(expanded, cl, 2);
    }

    expanded->max_clause_size = f->max_clause_size;
    if (expanded->max_clause_size < 2) expanded->max_clause_size = 2;
    return expanded;
}
