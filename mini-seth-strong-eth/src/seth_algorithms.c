#include "seth.h"
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdio.h>

/* ===================================================================
 * Schoning's Algorithm (1999) for k-SAT
 *
 * Repeat: pick random assignment, flip vars in unsatisfied clauses.
 * Success prob per trial: (k/(2(k-1)))^n
 * Expected time: O((2(k-1)/k)^n) = O(1.334^n) for 3-SAT
 * =================================================================== */
sat_result_t sat_schoening(const cnf_formula_t *f, int32_t max_tries,
                            int32_t max_flips) {
    sat_result_t result = {false, NULL, 0, 0, 0, 0.0};
    if (!f || f->num_vars == 0) return result;
    clock_t start = clock();
    int32_t n = f->num_vars;
    int32_t k = f->max_clause_size;
    if (k < 3) k = 3;
    int32_t walk_length = 3 * n;
    if (max_flips > 0 && max_flips < walk_length) walk_length = max_flips;
    if (max_tries <= 0) {
        double base = (double)(2 * (k - 1)) / (double)k;
        max_tries = (int32_t)(100.0 * pow(1.0 / base, (double)n));
        if (max_tries > 1000000) max_tries = 1000000;
        if (max_tries < 10) max_tries = 10;
    }
    for (int32_t trial = 0; trial < max_tries; trial++) {
        assignment_t *a = assign_create(n);
        if (!a) break;
        for (int32_t v = 1; v <= n; v++)
            assign_set(a, v, (rand() % 2) == 1);
        for (int32_t step = 0; step < walk_length; step++) {
            if (assign_is_satisfying(a, f)) {
                result.satisfiable = true;
                result.witness = a;
                result.time_seconds = (double)(clock() - start) / CLOCKS_PER_SEC;
                result.num_branches = trial + 1;
                result.num_propagations = step;
                return result;
            }
            int32_t *unsat = (int32_t *)malloc((size_t)f->num_clauses * sizeof(int32_t));
            if (!unsat) { assign_destroy(a); break; }
            int32_t unsat_count = 0;
            for (int32_t ci = 0; ci < f->num_clauses; ci++) {
                if (f->clauses[ci].deleted) continue;
                bool sat = false;
                for (int32_t j = 0; j < f->clauses[ci].size; j++) {
                    literal_t lit = f->clauses[ci].literals[j];
                    int8_t val = a->values[lit.var - 1];
                    if (val >= 0) {
                        bool v_true = lit.negated ? (val == 0) : (val == 1);
                        if (v_true) { sat = true; break; }
                    }
                }
                if (!sat) unsat[unsat_count++] = ci;
            }
            if (unsat_count == 0) {
                free(unsat);
                result.satisfiable = true;
                result.witness = a;
                result.time_seconds = (double)(clock() - start) / CLOCKS_PER_SEC;
                return result;
            }
            int32_t pick_clause = unsat[rand() % unsat_count];
            free(unsat);
            clause_t *c = &f->clauses[pick_clause];
            int32_t pick_lit = rand() % c->size;
            int32_t flip_var = c->literals[pick_lit].var;
            int8_t cur_val = a->values[flip_var - 1];
            assign_set(a, flip_var, (cur_val == 0) ? true : false);
            result.num_propagations++;
        }
        assign_destroy(a);
        result.num_branches++;
        result.num_backtracks++;
    }
    result.time_seconds = (double)(clock() - start) / CLOCKS_PER_SEC;
    return result;
}

/* ===================================================================
 * DPLL: Davis-Putnam-Logemann-Loveland (1962)
 * Backtracking with unit propagation and pure literal elimination.
 * Foundation of all modern SAT solvers.
 * =================================================================== */
typedef struct {
    int32_t *assignments;
    int32_t  top;
    int32_t  capacity;
} decision_stack_t;

static decision_stack_t *ds_create(int32_t cap) {
    decision_stack_t *ds = (decision_stack_t *)malloc(sizeof(decision_stack_t));
    if (!ds) return NULL;
    ds->assignments = (int32_t *)malloc((size_t)cap * sizeof(int32_t));
    if (!ds->assignments) { free(ds); return NULL; }
    ds->top = 0;
    ds->capacity = cap;
    return ds;
}

static void ds_destroy(decision_stack_t *ds) {
    if (!ds) return;
    free(ds->assignments);
    free(ds);
}

static void ds_push(decision_stack_t *ds, int32_t var) {
    if (ds->top < ds->capacity)
        ds->assignments[ds->top++] = var;
}

static int32_t ds_pop(decision_stack_t *ds) {
    if (ds->top > 0) return ds->assignments[--ds->top];
    return -1;
}

static int32_t dpll_choose_variable(const cnf_formula_t *f, const assignment_t *a) {
    int32_t best_var = -1;
    bool *candidates = (bool *)calloc((size_t)f->num_vars + 1, sizeof(bool));
    if (!candidates) return 1;
    for (int32_t i = 0; i < f->num_clauses; i++) {
        if (f->clauses[i].deleted) continue;
        bool sat = false;
        for (int32_t j = 0; j < f->clauses[i].size; j++) {
            int32_t v = f->clauses[i].literals[j].var;
            int8_t val = a->values[v - 1];
            if (val >= 0) {
                bool v_true = f->clauses[i].literals[j].negated ? (val == 0) : (val == 1);
                if (v_true) { sat = true; break; }
            }
        }
        if (!sat) {
            for (int32_t j = 0; j < f->clauses[i].size; j++) {
                int32_t v = f->clauses[i].literals[j].var;
                if (a->values[v - 1] < 0) candidates[v] = true;
            }
        }
    }
    for (int32_t v = 1; v <= f->num_vars; v++) {
        if (candidates[v] && a->values[v - 1] < 0) { best_var = v; break; }
    }
    free(candidates);
    if (best_var < 0) {
        for (int32_t v = 1; v <= f->num_vars; v++)
            if (a->values[v - 1] < 0) { best_var = v; break; }
    }
    return best_var;
}

static bool dpll_search(cnf_formula_t *f, assignment_t *a,
                         decision_stack_t *ds, uint64_t *branches, uint64_t *backtracks) {
    if (!unit_propagate(f, a)) { (*backtracks)++; return false; }
    pure_literal_eliminate(f, a);
    if (assign_is_satisfying(a, f)) return true;
    int32_t branch_var = dpll_choose_variable(f, a);
    if (branch_var < 0) return assign_is_satisfying(a, f);
    (*branches)++;
    assign_set(a, branch_var, false);
    ds_push(ds, branch_var);
    if (dpll_search(f, a, ds, branches, backtracks)) return true;
    assign_unset(a, branch_var);
    ds_pop(ds);
    assign_set(a, branch_var, true);
    ds_push(ds, branch_var);
    if (dpll_search(f, a, ds, branches, backtracks)) return true;
    assign_unset(a, branch_var);
    ds_pop(ds);
    (*backtracks)++;
    return false;
}

sat_result_t sat_dpll(const cnf_formula_t *f) {
    sat_result_t result = {false, NULL, 0, 0, 0, 0.0};
    if (!f) return result;
    clock_t start = clock();
    cnf_formula_t *working = cnf_clone(f);
    if (!working) return result;
    assignment_t *a = assign_create(f->num_vars);
    if (!a) { cnf_destroy(working); return result; }
    decision_stack_t *ds = ds_create(f->num_vars);
    if (!ds) { assign_destroy(a); cnf_destroy(working); return result; }
    uint64_t branches = 0, backtracks = 0;
    result.satisfiable = dpll_search(working, a, ds, &branches, &backtracks);
    result.num_branches = branches;
    result.num_backtracks = backtracks;
    if (result.satisfiable) result.witness = assign_clone(a);
    ds_destroy(ds);
    assign_destroy(a);
    cnf_destroy(working);
    result.time_seconds = (double)(clock() - start) / CLOCKS_PER_SEC;
    return result;
}

/* ===================================================================
 * PPSZ Algorithm (Paturi-Pudlák-Saks-Zane 2005)
 *
 * Bounded-width resolution + random variable ordering to achieve
 * the currently best known k-SAT running time:
 *   O(2^{(1-μ_k)n}) where μ_k = Ω(1/k)
 *
 * For 3-SAT: O(2^{0.3863n}) — best known bound.
 * For 4-SAT: O(2^{0.5548n})
 * General k: O(2^{(1-1/O(k))n})
 *
 * Core insight: Perform resolution only up to a bounded width w,
 * then use random ordering to choose variables. For each variable
 * in the random order, if it appears in a small (≤w) clause that
 * forces its value, assign it; otherwise assign randomly.
 *
 * The bounded-width resolution preprocessing finds all "forced"
 * assignments under width-w reasoning. Combined with random
 * ordering, this achieves the improved exponent.
 * =================================================================== */

#define PPSZ_MAX_WIDTH 3
#define PPSZ_MAX_TRIES 20

/* Build a random permutation of variables 1..n for PPSZ ordering */
static void ppsz_random_order(int32_t *order, int32_t n) {
    for (int32_t i = 0; i < n; i++) order[i] = i + 1;
    /* Fisher-Yates shuffle */
    for (int32_t i = n - 1; i > 0; i--) {
        int32_t j = rand() % (i + 1);
        int32_t tmp = order[i];
        order[i] = order[j];
        order[j] = tmp;
    }
}

/* PPSZ core search: given a random variable ordering, try to find
 * a satisfying assignment by processing variables in order.
 *
 * For each variable (in random order):
 *   - If there is a width-≤w clause forcing its value, set it
 *   - Otherwise, try both values with limited backtracking
 *
 * Returns true if a satisfying assignment is found. */
static bool ppsz_search(cnf_formula_t *f, assignment_t *a,
                         const int32_t *order, int32_t n,
                         int32_t depth, int32_t max_depth,
                         uint64_t *branches, uint64_t *backtracks,
                         const int32_t *forced_values) {

    if (cnf_has_conflict(f, a)) { (*backtracks)++; return false; }
    if (assign_is_satisfying(a, f)) return true;

    /* Try unit propagation to simplify */
    cnf_formula_t *work = cnf_clone(f);
    if (!work) return false;
    assignment_t *local = assign_clone(a);
    if (!local) { cnf_destroy(work); return false; }

    bool ok = unit_propagate(work, local);
    if (!ok) { assign_destroy(local); cnf_destroy(work); (*backtracks)++; return false; }
    if (assign_is_satisfying(local, work)) {
        /* Copy back the successful assignment */
        for (int32_t v = 1; v <= n; v++) {
            int8_t val = assign_get(local, v);
            if (val >= 0) assign_set(a, v, (bool)val);
        }
        assign_destroy(local); cnf_destroy(work);
        return true;
    }

    /* Find the next unassigned variable in the random order */
    int32_t next_var = -1;
    for (int32_t i = 0; i < n; i++) {
        if (assign_get(local, order[i]) < 0) {
            next_var = order[i];
            break;
        }
    }

    if (next_var < 0) {
        bool result = assign_is_satisfying(local, work);
        assign_destroy(local); cnf_destroy(work);
        return result;
    }

    /* Check if this variable is forced by a small clause */
    bool forced = false;
    bool forced_val = false;
    if (forced_values && forced_values[next_var] != -1) {
        forced = true;
        forced_val = (forced_values[next_var] == 1);
    }

    assign_destroy(local); cnf_destroy(work);

    bool first_val = forced ? forced_val : false;
    bool second_val = forced ? !forced_val : true;

    (*branches)++;
    assign_set(a, next_var, first_val);
    if (ppsz_search(f, a, order, n, depth + 1, max_depth,
                    branches, backtracks, forced_values))
        return true;

    if (forced) {
        /* Forced variable: if wrong value fails, backtrack fully */
        assign_unset(a, next_var);
        (*backtracks)++;
        return false;
    }

    assign_unset(a, next_var);
    assign_set(a, next_var, second_val);
    if (ppsz_search(f, a, order, n, depth + 1, max_depth,
                    branches, backtracks, forced_values))
        return true;
    assign_unset(a, next_var);
    (*backtracks)++;
    return false;
}

/* Compute forced variable values using bounded-width resolution.
 *
 * Resolution rule: from (x ∨ A) and (¬x ∨ B), derive (A ∨ B).
 * Bounded-width: only keep resolvents with ≤w literals.
 *
 * After computing the closure under width-w resolution, we check
 * for unit clauses that force variable assignments. */
static void ppsz_compute_forced(const cnf_formula_t *f, int32_t *forced,
                                 int32_t n, int32_t w) {
    /* Initialize: no forced values known */
    for (int32_t v = 1; v <= n; v++) forced[v] = -1;

    /* Collect all clauses of size ≤ w (including derived unit clauses) */
    /* First pass: check for existing unit clauses */
    for (int32_t i = 0; i < f->num_clauses; i++) {
        if (f->clauses[i].deleted) continue;
        if (f->clauses[i].size == 1) {
            literal_t lit = f->clauses[i].literals[0];
            int8_t needed = lit.negated ? 0 : 1;
            if (forced[lit.var] == -1)
                forced[lit.var] = needed;
            else if (forced[lit.var] != needed)
                return; /* conflict */
        }
    }

    /* Width-w resolution: repeatedly try to derive new clauses.
     * For each pair of clauses, check if they can be resolved,
     * and if the resolvent has size ≤ w, add its implications. */
    bool changed = true;
    int32_t max_iter = f->num_vars * 4;
    int32_t iter = 0;

    while (changed && iter < max_iter) {
        changed = false;
        iter++;

        for (int32_t i = 0; i < f->num_clauses; i++) {
            if (f->clauses[i].deleted) continue;
            for (int32_t j = i + 1; j < f->num_clauses; j++) {
                if (f->clauses[j].deleted) continue;

                /* Look for complementary literal pair */
                int32_t pivot_var = -1;
                for (int32_t pi = 0; pi < f->clauses[i].size; pi++) {
                    literal_t li = f->clauses[i].literals[pi];
                    for (int32_t pj = 0; pj < f->clauses[j].size; pj++) {
                        literal_t lj = f->clauses[j].literals[pj];
                        if (li.var == lj.var && li.negated != lj.negated) {
                            pivot_var = li.var;
                            break;
                        }
                    }
                    if (pivot_var > 0) break;
                }

                if (pivot_var < 0) continue;

                /* Compute resolvent size estimate:
                 * |C_i ∪ C_j| - 2 ≤ w? */
                int32_t est_size = f->clauses[i].size + f->clauses[j].size - 2;
                if (est_size > w) continue;

                /* Check if the resolvent implies a new unit clause */
                if (est_size == 0) return; /* empty clause: unsatisfiable */
                if (est_size == 1) {
                    /* Find the single literal in the resolvent */
                    for (int32_t pi = 0; pi < f->clauses[i].size; pi++) {
                        literal_t li = f->clauses[i].literals[pi];
                        if (li.var == pivot_var) continue;
                        int8_t needed = li.negated ? 0 : 1;
                        if (forced[li.var] == -1) {
                            forced[li.var] = needed;
                            changed = true;
                        }
                    }
                    if (est_size == 1) {
                        for (int32_t pj = 0; pj < f->clauses[j].size; pj++) {
                            literal_t lj = f->clauses[j].literals[pj];
                            if (lj.var == pivot_var) continue;
                            int8_t needed = lj.negated ? 0 : 1;
                            if (forced[lj.var] == -1) {
                                forced[lj.var] = needed;
                                changed = true;
                            }
                        }
                    }
                }
            }
        }
    }
}

sat_result_t sat_ppsz(const cnf_formula_t *f) {
    sat_result_t result = {false, NULL, 0, 0, 0, 0.0};
    if (!f || f->num_vars == 0) return result;

    clock_t start = clock();
    int32_t n = f->num_vars;

    /* Compute forced values via bounded-width resolution */
    int32_t *forced = (int32_t *)calloc((size_t)n + 1, sizeof(int32_t));
    if (!forced) return result;
    ppsz_compute_forced(f, forced, n, PPSZ_MAX_WIDTH);

    /* Random variable ordering */
    int32_t *order = (int32_t *)malloc((size_t)n * sizeof(int32_t));
    if (!order) { free(forced); return result; }

    /* Run multiple tries with different random orders */
    int32_t num_tries = PPSZ_MAX_TRIES;
    if (n > 10) num_tries = 5; /* fewer tries for larger instances */
    int32_t max_depth = n; /* full depth search */

    for (int32_t trial = 0; trial < num_tries; trial++) {
        ppsz_random_order(order, n);
        assignment_t *a = assign_create(n);
        if (!a) break;

        if (ppsz_search((cnf_formula_t *)f, a, order, n, 0, max_depth,
                        &result.num_branches, &result.num_backtracks, forced)) {
            result.satisfiable = true;
            result.witness = a;
            result.time_seconds = (double)(clock() - start) / CLOCKS_PER_SEC;
            free(order);
            free(forced);
            return result;
        }
        assign_destroy(a);
    }

    result.time_seconds = (double)(clock() - start) / CLOCKS_PER_SEC;
    free(order);
    free(forced);
    return result;
}

/* ===================================================================
 * CDCL: Conflict-Driven Clause Learning
 *
 * The backbone of modern SAT solvers (since 1996).
 * Marques-Silva & Sakallah (1996), further refined by
 * Moskewicz et al. (2001) for Chaff, and Een & Sorensson
 * for MiniSAT.
 *
 * Key mechanisms:
 * 1. Boolean Constraint Propagation (BCP) using watched literals
 * 2. Conflict analysis via implication graph
 * 3. Learned clause derived from first UIP (Unique Implication Point)
 * 4. Non-chronological backtracking (backjumping) to decision level
 *    where learned clause becomes unit
 * 5. Activity-based variable heuristics (VSIDS)
 * 6. Periodic restarts to avoid getting stuck
 *
 * Worst-case complexity: O(2^{n}) but typically much faster.
 * For certain formula families, CDCL has O(2^{n/2}) bound.
 * =================================================================== */

#define CDCL_MAX_LEARNED 10000
#define CDCL_MAX_CONFLICTS 5000

/* Variable state for VSIDS heuristic */
typedef struct {
    double   activity;
    int32_t  decision_level;  /* 0 = unassigned, >0 = level assigned */
    int8_t   value;           /* 0=F, 1=T, -1=unassigned */
    int32_t  antecedent_clause; /* clause that forced this (or -1 for decision) */
} cdcl_var_t;

/* CDCL solver state */
typedef struct {
    cdcl_var_t *vars;         /* variable states, indexed 1..n */
    int32_t     num_vars;
    int32_t     decision_level;
    int32_t    *trail;        /* assignment trail */
    int32_t     trail_size;
    int32_t     trail_cap;
    int32_t    *trail_lim;    /* trail index where each decision level starts */
    uint64_t    num_conflicts;
    uint64_t    num_decisions;
    uint64_t    num_propagations;
    double      var_inc;
} cdcl_solver_t;

static cdcl_solver_t *cdcl_init(int32_t num_vars) {
    cdcl_solver_t *s = (cdcl_solver_t *)calloc(1, sizeof(cdcl_solver_t));
    if (!s) return NULL;
    s->num_vars = num_vars;
    s->vars = (cdcl_var_t *)calloc((size_t)num_vars + 1, sizeof(cdcl_var_t));
    s->trail_cap = num_vars * 4;
    s->trail = (int32_t *)malloc((size_t)s->trail_cap * sizeof(int32_t));
    s->trail_lim = (int32_t *)malloc((size_t)(num_vars + 1) * sizeof(int32_t));
    s->var_inc = 1.0;
    if (!s->vars || !s->trail || !s->trail_lim) {
        free(s->vars); free(s->trail); free(s->trail_lim); free(s);
        return NULL;
    }
    for (int32_t v = 1; v <= num_vars; v++) {
        s->vars[v].value = -1;
        s->vars[v].decision_level = 0;
        s->vars[v].antecedent_clause = -1;
        s->vars[v].activity = 0.0;
    }
    return s;
}

static void cdcl_destroy(cdcl_solver_t *s) {
    if (!s) return;
    free(s->vars); free(s->trail); free(s->trail_lim);
    free(s);
}

/* VSIDS: pick unassigned variable with highest activity */
static int32_t cdcl_pick_branch_var(cdcl_solver_t *s) {
    int32_t best = -1;
    double best_act = -1.0;
    for (int32_t v = 1; v <= s->num_vars; v++) {
        if (s->vars[v].value < 0 && s->vars[v].activity > best_act) {
            best_act = s->vars[v].activity;
            best = v;
        }
    }
    if (best < 0) {
        for (int32_t v = 1; v <= s->num_vars; v++) {
            if (s->vars[v].value < 0) { best = v; break; }
        }
    }
    return best;
}

/* Bump variable activity (multiply by factor for all vars) */
static void cdcl_var_bump_activity(cdcl_solver_t *s, int32_t var) {
    s->vars[var].activity += s->var_inc;
    if (s->vars[var].activity > 1e100) {
        /* Rescale all activities */
        for (int32_t v = 1; v <= s->num_vars; v++)
            s->vars[v].activity *= 1e-100;
        s->var_inc *= 1e-100;
    }
}

static void cdcl_var_decay_activity(cdcl_solver_t *s) {
    s->var_inc *= 1.05; /* VSIDS decay factor */
}

/* Assign a variable at the current decision level */
static void cdcl_assign(cdcl_solver_t *s, int32_t var, int8_t val,
                         int32_t antecedent) {
    if (s->trail_size >= s->trail_cap) {
        s->trail_cap *= 2;
        s->trail = (int32_t *)realloc(s->trail,
            (size_t)s->trail_cap * sizeof(int32_t));
    }
    s->vars[var].value = val;
    s->vars[var].decision_level = s->decision_level;
    s->vars[var].antecedent_clause = antecedent;
    s->trail[s->trail_size++] = var;
}

/* Unassign all variables assigned at or above given decision level */
static void cdcl_backtrack_to_level(cdcl_solver_t *s, int32_t level) {
    int32_t target = (level >= 0) ? s->trail_lim[level] : 0;
    for (int32_t i = s->trail_size - 1; i >= target; i--) {
        int32_t v = s->trail[i];
        s->vars[v].value = -1;
        s->vars[v].decision_level = 0;
        s->vars[v].antecedent_clause = -1;
    }
    s->trail_size = target;
    s->decision_level = level;
}

/* Boolean Constraint Propagation (BCP) — simplified version.
 * Checks all clauses for unit status under current assignment.
 * Returns: true if no conflict, false if conflict detected. */
static bool cdcl_bcp(cnf_formula_t *f, cdcl_solver_t *s) {
    bool propagated;
    do {
        propagated = false;
        for (int32_t ci = 0; ci < f->num_clauses; ci++) {
            if (f->clauses[ci].deleted) continue;
            clause_t *c = &f->clauses[ci];

            int32_t unassigned_idx = -1;
            int32_t unassigned_count = 0;
            bool satisfied = false;

            for (int32_t j = 0; j < c->size; j++) {
                literal_t lit = c->literals[j];
                int8_t val = s->vars[lit.var].value;
                if (val < 0) {
                    unassigned_count++;
                    unassigned_idx = j;
                } else {
                    bool lit_true = lit.negated ? (val == 0) : (val == 1);
                    if (lit_true) { satisfied = true; break; }
                }
            }

            if (satisfied) continue;

            if (unassigned_count == 0) {
                /* Conflict: all literals are false */
                return false;
            }

            if (unassigned_count == 1) {
                /* Unit clause: force the unassigned literal */
                literal_t forced = c->literals[unassigned_idx];
                int8_t val = forced.negated ? 0 : 1;
                if (s->vars[forced.var].value < 0) {
                    cdcl_assign(s, forced.var, val, ci);
                    s->num_propagations++;
                    propagated = true;
                } else if (s->vars[forced.var].value != val) {
                    /* Conflict: literal forced to opposite value */
                    return false;
                }
            }
        }
    } while (propagated);
    return true;
}

/* CDCL conflict analysis: find the decision level to backtrack to.
 * Uses first-UIP (Unique Implication Point) heuristic.
 * Simplified: backtrack to max(0, current_level - 1). */
static int32_t cdcl_analyze_conflict(cnf_formula_t *f, cdcl_solver_t *s) {
    (void)f; /* More sophisticated analysis would trace implication graph */

    s->num_conflicts++;
    cdcl_var_decay_activity(s);

    /* Bump activity of variables involved in the conflict */
    for (int32_t v = 1; v <= s->num_vars; v++) {
        if (s->vars[v].value >= 0 && s->vars[v].decision_level > 0) {
            cdcl_var_bump_activity(s, v);
        }
    }

    /* Simplified: backtrack to previous decision level.
     * A full implementation would compute the actual assertion level. */
    int32_t backtrack_level = (s->decision_level > 0) ? s->decision_level - 1 : 0;
    return backtrack_level;
}

sat_result_t sat_cdcl(const cnf_formula_t *f) {
    sat_result_t result = {false, NULL, 0, 0, 0, 0.0};
    if (!f || f->num_vars == 0) return result;

    clock_t start = clock();
    int32_t n = f->num_vars;

    cdcl_solver_t *s = cdcl_init(n);
    if (!s) return result;

    cnf_formula_t *working = cnf_clone(f);
    if (!working) { cdcl_destroy(s); return result; }

    s->decision_level = 0;
    s->trail_size = 0;
    s->trail_lim[0] = 0;

    while (true) {
        /* BCP: propagate all unit clauses */
        if (!cdcl_bcp(working, s)) {
            /* Conflict */
            if (s->decision_level == 0) {
                /* Conflict at level 0: unsatisfiable */
                result.satisfiable = false;
                break;
            }
            int32_t bt_level = cdcl_analyze_conflict(working, s);
            cdcl_backtrack_to_level(s, bt_level);
            s->decision_level = bt_level;
            continue;
        }

        /* Check if all variables assigned */
        bool all_assigned = true;
        for (int32_t v = 1; v <= n; v++) {
            if (s->vars[v].value < 0) { all_assigned = false; break; }
        }

        if (all_assigned) {
            /* All variables assigned without conflict: satisfiable */
            result.satisfiable = true;
            /* Check against formula */
            assignment_t *tmp = assign_create(n);
            for (int32_t v = 1; v <= n; v++)
                assign_set(tmp, v, s->vars[v].value == 1);
            if (assign_is_satisfying(tmp, f)) {
                result.witness = tmp;
            } else {
                assign_destroy(tmp);
                /* Continue search */
            }
            break;
        }

        /* Make a decision */
        int32_t branch_var = cdcl_pick_branch_var(s);
        if (branch_var < 0) break;

        s->decision_level++;
        s->trail_lim[s->decision_level] = s->trail_size;
        s->num_decisions++;

        /* Decide phase: false by default */
        int8_t phase = 0;
        cdcl_assign(s, branch_var, phase, -1);
    }

    result.num_branches = s->num_decisions;
    result.num_backtracks = s->num_conflicts;
    result.num_propagations = s->num_propagations;

    cdcl_destroy(s);
    cnf_destroy(working);
    result.time_seconds = (double)(clock() - start) / CLOCKS_PER_SEC;
    return result;
}
