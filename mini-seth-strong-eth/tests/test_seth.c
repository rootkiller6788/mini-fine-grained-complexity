/* test_seth.c -- Tests for mini-seth-strong-eth submodule */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "../include/seth.h"
#include "../include/cnf_generator.h"
#include "../include/reduction.h"
#include "../include/exponential_complexity.h"

static int tests_run = 0;
static int tests_passed = 0;
#define TEST(name) do { tests_run++; printf("  TEST %s ... ", name); } while(0)
#define PASS() do { tests_passed++; printf("PASS\n"); } while(0)

static void test_cnf_create_destroy(void) {
    TEST("cnf_create_destroy");
    cnf_formula_t *f = cnf_create(5, 10);
    assert(f != NULL);
    assert(f->num_vars == 5);
    cnf_destroy(f);
    PASS();
}

static void test_add_clause(void) {
    TEST("add_clause");
    cnf_formula_t *f = cnf_create(3, 10);
    int32_t c1[] = {1, -2, 3};
    cnf_add_clause(f, c1, 3);
    assert(f->num_clauses == 1);
    assert(f->clauses[0].size == 3);
    cnf_destroy(f);
    PASS();
}

static void test_assignment(void) {
    TEST("assignment_ops");
    assignment_t *a = assign_create(5);
    assert(a != NULL);
    assign_set(a, 1, true);
    assert(assign_get(a, 1) == 1);
    assign_unset(a, 1);
    assert(assign_get(a, 1) == -1);
    assign_destroy(a);
    PASS();
}

static void test_satisfiability_check(void) {
    TEST("satisfiability_check");
    cnf_formula_t *f = cnf_create(2, 10);
    int32_t c1[] = {1, 2};
    int32_t c2[] = {-1, -2};
    cnf_add_clause(f, c1, 2);
    cnf_add_clause(f, c2, 2);
    assignment_t *a = assign_create(2);
    assign_set(a, 1, true);
    assign_set(a, 2, false);
    assert(assign_is_satisfying(a, f));
    assign_destroy(a);
    cnf_destroy(f);
    PASS();
}

static void test_seth_limit(void) {
    TEST("seth_limit");
    double s3 = seth_limit_s_k(3);
    assert(s3 > 0.3 && s3 < 0.5);
    PASS();
}

static void test_unit_propagation(void) {
    TEST("unit_propagation");
    cnf_formula_t *f = cnf_create(2, 10);
    int32_t c1[] = {1};
    int32_t c2[] = {-1, 2};
    cnf_add_clause(f, c1, 1);
    cnf_add_clause(f, c2, 2);
    assignment_t *a = assign_create(2);
    bool ok = unit_propagate(f, a);
    assert(ok);
    assert(assign_get(a, 1) == 1);
    assert(assign_get(a, 2) == 1);
    assign_destroy(a);
    cnf_destroy(f);
    PASS();
}

static void test_unit_conflict(void) {
    TEST("unit_conflict");
    cnf_formula_t *f = cnf_create(1, 10);
    int32_t c1[] = {1};
    int32_t c2[] = {-1};
    cnf_add_clause(f, c1, 1);
    cnf_add_clause(f, c2, 1);
    assignment_t *a = assign_create(1);
    bool ok = unit_propagate(f, a);
    assert(!ok);
    assign_destroy(a);
    cnf_destroy(f);
    PASS();
}

static void test_brute_force_sat(void) {
    TEST("brute_force_sat");
    cnf_formula_t *f = cnf_create(3, 10);
    int32_t c1[] = {1, 2};
    int32_t c2[] = {-1, 3};
    cnf_add_clause(f, c1, 2);
    cnf_add_clause(f, c2, 2);
    sat_result_t r = sat_brute_force(f);
    assert(r.satisfiable);
    if (r.witness) assign_destroy(r.witness);
    cnf_destroy(f);
    PASS();
}

static void test_dpll_sat(void) {
    TEST("dpll_sat");
    cnf_formula_t *f = cnf_create(3, 10);
    int32_t c1[] = {1, 2, 3};
    int32_t c2[] = {-1, -2};
    cnf_add_clause(f, c1, 3);
    cnf_add_clause(f, c2, 2);
    sat_result_t r = sat_dpll(f);
    assert(r.satisfiable);
    if (r.witness) assign_destroy(r.witness);
    cnf_destroy(f);
    PASS();
}

static void test_random_ksat(void) {
    TEST("random_ksat");
    cnf_formula_t *f = generate_random_ksat(5, 3, 4.0);
    assert(f != NULL);
    cnf_destroy(f);
    PASS();
}

static void test_pigeonhole(void) {
    TEST("pigeonhole_unsat");
    cnf_formula_t *f = generate_pigeonhole(2);
    assert(f != NULL);
    sat_result_t r = sat_brute_force(f);
    assert(!r.satisfiable);
    cnf_destroy(f);
    PASS();
}

static void test_sparsification(void) {
    TEST("sparsification");
    cnf_formula_t *f = generate_random_ksat(5, 3, 5.0);
    assert(f != NULL);
    sparsification_result_t *sr = sparsify(f, 0.5);
    assert(sr != NULL);
    assert(sr->num_formulas > 0);
    sparsification_free(sr);
    cnf_destroy(f);
    PASS();
}

static void test_reduction_graph(void) {
    TEST("reduction_graph");
    reduction_graph_t *g = build_standard_reduction_graph();
    assert(g != NULL);
    assert(g->num_problems >= 5);
    reduction_graph_destroy(g);
    PASS();
}

int main(void) {
    printf("\n=== mini-seth-strong-eth Test Suite ===\n\n");
    test_cnf_create_destroy();
    test_add_clause();
    test_assignment();
    test_satisfiability_check();
    test_seth_limit();
    test_unit_propagation();
    test_unit_conflict();
    test_brute_force_sat();
    test_dpll_sat();
    test_random_ksat();
    test_pigeonhole();
    test_sparsification();
    test_reduction_graph();
    printf("\n=== Results: %d/%d tests passed ===\n\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}
