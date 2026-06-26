/* ============================================================================
 * seth.c ? Core SETH/ETH implementation
 *
 * Implements the Strong Exponential Time Hypothesis (SETH) and
 * Exponential Time Hypothesis (ETH) formal framework.
 *
 * Knowledge coverage:
 *   L1: struct definitions for CNF, assignment, literal
 *   L2: SETH/ETH violation checking, fine-grained reduction concepts
 *   L3: CNF formula mathematical structure
 *   L4: Sparsification Lemma, ETH?SETH relationship
 *   L5: Brute-force SAT, unit propagation, pure literal elimination
 *
 * References:
 *   - Impagliazzo, Paturi, Zane (2001) "Which Problems Have Strongly
 *     Exponential Complexity?" JCSS 63(4):512-530
 *   - Impagliazzo, Paturi (1999) "Complexity of k-SAT"
 *   - Sch?ning (1999) "A Probabilistic Algorithm for k-SAT" FOCS
 * ============================================================================ */

#include "seth.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>

/* ============================================================================
 * L1: Core CNF/Assignment/Literal Operations
 * ============================================================================ */

/* Convert a signed integer literal to literal_t */
static literal_t literal_from_int(int32_t signed_lit) {
    literal_t l;
    l.var = (signed_lit > 0) ? signed_lit : -signed_lit;
    l.negated = (signed_lit < 0);
    return l;
}

/* Evaluate a literal under a given assignment */
static bool literal_eval(literal_t l, const assignment_t *a) {
    if (l.var < 1 || l.var > a->num_vars) return false;
    int8_t val = a->values[l.var - 1];
    if (val < 0) return false;
    return l.negated ? (val == 0) : (val == 1);
}

/* Check if literal is falsified under assignment */
static bool literal_falsified(literal_t l, const assignment_t *a) {
    if (l.var < 1 || l.var > a->num_vars) return false;
    int8_t val = a->values[l.var - 1];
    if (val < 0) return false;
    return l.negated ? (val == 1) : (val == 0);
}

/* ---- CNF Formula Operations ---- */

cnf_formula_t *cnf_create(int32_t num_vars, int32_t capacity) {
    cnf_formula_t *f = (cnf_formula_t *)malloc(sizeof(cnf_formula_t));
    if (!f) return NULL;
    f->num_vars = num_vars;
    f->num_clauses = 0;
    f->max_clause_size = 0;
    f->capacity = (capacity > 0) ? capacity : 64;
    f->clauses = (clause_t *)calloc((size_t)f->capacity, sizeof(clause_t));
    if (!f->clauses) { free(f); return NULL; }
    return f;
}

void cnf_destroy(cnf_formula_t *f) {
    if (!f) return;
    for (int32_t i = 0; i < f->num_clauses; i++) {
        free(f->clauses[i].literals);
    }
    free(f->clauses);
    free(f);
}

clause_t cnf_add_clause(cnf_formula_t *f, int32_t *lits, int32_t n) {
    if (f->num_clauses >= f->capacity) {
        int32_t new_cap = f->capacity * 2;
        clause_t *new_clauses = (clause_t *)realloc(f->clauses,
            (size_t)new_cap * sizeof(clause_t));
        if (!new_clauses) {
            clause_t empty = {NULL, 0, false};
            return empty;
        }
        f->clauses = new_clauses;
        for (int32_t i = f->capacity; i < new_cap; i++) {
            f->clauses[i].literals = NULL;
            f->clauses[i].size = 0;
            f->clauses[i].deleted = false;
        }
        f->capacity = new_cap;
    }
    clause_t *c = &f->clauses[f->num_clauses];
    c->size = n;
    c->deleted = false;
    c->literals = (literal_t *)malloc((size_t)n * sizeof(literal_t));
    if (!c->literals) {
        clause_t empty = {NULL, 0, false};
        return empty;
    }
    for (int32_t i = 0; i < n; i++) {
        c->literals[i] = literal_from_int(lits[i]);
    }
    if (n > f->max_clause_size) f->max_clause_size = n;
    f->num_clauses++;
    return *c;
}

void cnf_add_clause_raw(cnf_formula_t *f, literal_t *literals, int32_t n) {
    if (f->num_clauses >= f->capacity) {
        int32_t new_cap = f->capacity * 2;
        f->clauses = (clause_t *)realloc(f->clauses,
            (size_t)new_cap * sizeof(clause_t));
        for (int32_t i = f->capacity; i < new_cap; i++) {
            f->clauses[i].literals = NULL;
            f->clauses[i].size = 0;
            f->clauses[i].deleted = false;
        }
        f->capacity = new_cap;
    }
    clause_t *c = &f->clauses[f->num_clauses];
    c->size = n;
    c->deleted = false;
    c->literals = (literal_t *)malloc((size_t)n * sizeof(literal_t));
    if (!c->literals) return;
    memcpy(c->literals, literals, (size_t)n * sizeof(literal_t));
    if (n > f->max_clause_size) f->max_clause_size = n;
    f->num_clauses++;
}

cnf_formula_t *cnf_clone(const cnf_formula_t *f) {
    cnf_formula_t *clone = cnf_create(f->num_vars, f->capacity);
    if (!clone) return NULL;
    for (int32_t i = 0; i < f->num_clauses; i++) {
        if (!f->clauses[i].deleted) {
            cnf_add_clause_raw(clone, f->clauses[i].literals,
                               f->clauses[i].size);
        }
    }
    clone->max_clause_size = f->max_clause_size;
    return clone;
}

int32_t cnf_count_clauses(const cnf_formula_t *f) {
    int32_t count = 0;
    for (int32_t i = 0; i < f->num_clauses; i++) {
        if (!f->clauses[i].deleted) count++;
    }
    return count;
}

/* ---- Assignment Operations ---- */

assignment_t *assign_create(int32_t num_vars) {
    assignment_t *a = (assignment_t *)malloc(sizeof(assignment_t));
    if (!a) return NULL;
    a->num_vars = num_vars;
    a->num_assigned = 0;
    a->values = (int8_t *)malloc((size_t)num_vars * sizeof(int8_t));
    if (!a->values) { free(a); return NULL; }
    for (int32_t i = 0; i < num_vars; i++) {
        a->values[i] = -1;
    }
    return a;
}

void assign_destroy(assignment_t *a) {
    if (!a) return;
    free(a->values);
    free(a);
}

void assign_set(assignment_t *a, int32_t var, bool value) {
    if (var < 1 || var > a->num_vars) return;
    if (a->values[var - 1] < 0) a->num_assigned++;
    a->values[var - 1] = value ? 1 : 0;
}

void assign_unset(assignment_t *a, int32_t var) {
    if (var < 1 || var > a->num_vars) return;
    if (a->values[var - 1] >= 0) a->num_assigned--;
    a->values[var - 1] = -1;
}

int8_t assign_get(const assignment_t *a, int32_t var) {
    if (var < 1 || var > a->num_vars) return -1;
    return a->values[var - 1];
}

assignment_t *assign_clone(const assignment_t *a) {
    assignment_t *clone = assign_create(a->num_vars);
    if (!clone) return NULL;
    memcpy(clone->values, a->values, (size_t)a->num_vars * sizeof(int8_t));
    clone->num_assigned = a->num_assigned;
    return clone;
}

bool assign_is_satisfying(const assignment_t *a, const cnf_formula_t *f) {
    for (int32_t i = 0; i < f->num_clauses; i++) {
        if (f->clauses[i].deleted) continue;
        bool clause_sat = false;
        for (int32_t j = 0; j < f->clauses[i].size; j++) {
            if (literal_eval(f->clauses[i].literals[j], a)) {
                clause_sat = true;
                break;
            }
        }
        if (!clause_sat) return false;
    }
    return true;
}

bool cnf_is_satisfied(const cnf_formula_t *f, const assignment_t *a) {
    return assign_is_satisfying(a, f);
}

bool cnf_has_conflict(const cnf_formula_t *f, const assignment_t *a) {
    for (int32_t i = 0; i < f->num_clauses; i++) {
        if (f->clauses[i].deleted) continue;
        bool all_false = true;
        for (int32_t j = 0; j < f->clauses[i].size; j++) {
            if (!literal_falsified(f->clauses[i].literals[j], a)) {
                all_false = false;
                break;
            }
        }
        if (all_false) return true;
    }
    return false;
}

/* ============================================================================
 * L2: ETH/SETH Formal Statements
 * ============================================================================ */

bool seth_violated(double achieved_exponent, int32_t k) {
    if (k < 3) return false;
    double limit = seth_limit_s_k(k);
    return achieved_exponent < limit - 1e-9;
}

double seth_limit_s_k(int32_t k) {
    static const double s_values[] = {
        0.0, 0.0, 0.0,
        0.3863, 0.5548, 0.6500, 0.7162, 0.7627,
        0.7977, 0.8245, 0.8453
    };
    if (k <= 10) return s_values[k];
    return 1.0 - (1.5 / (double)k);
}

bool eth_violated(double achieved_exponent, int32_t n) {
    (void)n;
    return achieved_exponent < 0.001;
}

/* ============================================================================
 * Unit Propagation and Simplification (L5)
 * ============================================================================ */

bool unit_propagate(cnf_formula_t *f, assignment_t *a) {
    bool changed = true;
    int32_t iterations = 0;
    int32_t max_iterations = f->num_vars * f->num_clauses + 100;
    while (changed && iterations < max_iterations) {
        changed = false;
        iterations++;
        for (int32_t i = 0; i < f->num_clauses; i++) {
            if (f->clauses[i].deleted) continue;
            clause_t *c = &f->clauses[i];
            int32_t unassigned_count = 0;
            int32_t unassigned_idx = -1;
            bool any_satisfied = false;
            for (int32_t j = 0; j < c->size; j++) {
                literal_t lit = c->literals[j];
                int8_t val = a->values[lit.var - 1];
                if (val < 0) {
                    unassigned_count++;
                    unassigned_idx = j;
                } else {
                    bool lit_true = lit.negated ? (val == 0) : (val == 1);
                    if (lit_true) { any_satisfied = true; break; }
                }
            }
            if (any_satisfied) continue;
            if (unassigned_count == 0) return false;
            if (unassigned_count == 1) {
                literal_t forced = c->literals[unassigned_idx];
                assign_set(a, forced.var, !forced.negated);
                changed = true;
            }
        }
    }
    return true;
}

void pure_literal_eliminate(cnf_formula_t *f, assignment_t *a) {
    int32_t n = f->num_vars;
    int32_t *pos_count = (int32_t *)calloc((size_t)n + 1, sizeof(int32_t));
    int32_t *neg_count = (int32_t *)calloc((size_t)n + 1, sizeof(int32_t));
    if (!pos_count || !neg_count) {
        free(pos_count); free(neg_count);
        return;
    }
    for (int32_t i = 0; i < f->num_clauses; i++) {
        if (f->clauses[i].deleted) continue;
        for (int32_t j = 0; j < f->clauses[i].size; j++) {
            literal_t lit = f->clauses[i].literals[j];
            if (lit.negated) neg_count[lit.var]++;
            else pos_count[lit.var]++;
        }
    }
    for (int32_t v = 1; v <= n; v++) {
        if (a->values[v - 1] >= 0) continue;
        if (pos_count[v] > 0 && neg_count[v] == 0)
            assign_set(a, v, true);
        else if (neg_count[v] > 0 && pos_count[v] == 0)
            assign_set(a, v, false);
    }
    free(pos_count);
    free(neg_count);
}

/* ============================================================================
 * L5: Brute-force SAT solver (baseline for SETH lower bounds)
 * ============================================================================ */

static bool bf_recursive_inner(cnf_formula_t *f, assignment_t *a,
                                int32_t var, uint64_t *branches) {
    if (var > f->num_vars)
        return assign_is_satisfying(a, f);
    if (cnf_has_conflict(f, a)) return false;
    (*branches)++;
    assign_set(a, var, false);
    if (bf_recursive_inner(f, a, var + 1, branches)) return true;
    assign_unset(a, var);
    assign_set(a, var, true);
    if (bf_recursive_inner(f, a, var + 1, branches)) return true;
    assign_unset(a, var);
    return false;
}

sat_result_t sat_brute_force(const cnf_formula_t *f) {
    sat_result_t result = {false, NULL, 0, 0, 0, 0.0};
    if (!f) return result;
    clock_t start = clock();
    assignment_t *a = assign_create(f->num_vars);
    if (!a) return result;
    uint64_t branches = 0;
    result.satisfiable = bf_recursive_inner((cnf_formula_t *)f, a, 1, &branches);
    result.num_branches = branches;
    if (result.satisfiable)
        result.witness = assign_clone(a);
    assign_destroy(a);
    result.time_seconds = (double)(clock() - start) / CLOCKS_PER_SEC;
    return result;
}

/* PPSZ and CDCL are implemented in seth_algorithms.c.
 * The declarations here serve as the API entry points. */
