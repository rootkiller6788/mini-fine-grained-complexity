/* test_main.c -- Comprehensive tests for mini-conditional-lower-bounds */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "condlb.h"
#include "seth_bound.h"
#include "ov_reduction.h"
#include "threesum_apsp.h"
#include "string_hardness.h"
#include "reduction_web.h"

static int tests_run = 0;
static int tests_passed = 0;

#define TEST(name) do { tests_run++; printf("  TEST: %s ... ", name); } while(0)
#define PASS() do { printf("PASS\n"); tests_passed++; } while(0)
#define FAIL(msg) do { printf("FAIL: %s\n", msg); } while(0)

/* L1: Definitions tests */
static void test_hypothesis_db(void) {
    TEST("Hypothesis database create/free");
    HypothesisDatabase* db = hyp_db_create(10, 20);
    assert(db != NULL);
    assert(db->n_hypotheses == 0);
    assert(db->n_bounds == 0);
    PASS();

    TEST("Add hypothesis");
    int id = hyp_db_add_hypothesis(db, "SETH", "k-SAT", 2.0, 1e-9,
        HYPOTHESIS_UNREFUTED, 2001, "Impagliazzo & Paturi");
    assert(id == 0);
    assert(db->n_hypotheses == 1);
    PASS();

    TEST("Find hypothesis by name");
    int found = hyp_db_find_hypothesis(db, "SETH");
    assert(found == 0);
    assert(hyp_db_find_hypothesis(db, "NONEXISTENT") == -1);
    PASS();

    TEST("Add implication");
    int id2 = hyp_db_add_hypothesis(db, "ETH", "3-SAT", 1.0, 0.0,
        HYPOTHESIS_UNREFUTED, 1999, "Impagliazzo & Paturi");
    hyp_db_add_implication(db, id, id2);
    assert(db->hypotheses[id].n_implies == 1);
    assert(db->hypotheses[id2].n_depends_on == 1);
    PASS();

    TEST("Add conditional lower bound");
    int bid = hyp_db_add_bound(db, id, "Edit Distance", 2.0,
        REDUCTION_FINE_GRAINED, 1, "Backurs & Indyk (STOC 2015)");
    assert(bid == 0);
    assert(db->n_bounds == 1);
    PASS();

    hyp_db_free(db);
}

/* L3: Reduction web tests */
static void test_reduction_web(void) {
    TEST("Reduction web create/free");
    ReductionWeb* rw = rw_create(20);
    assert(rw != NULL);
    assert(rw->n_problems == 0);
    PASS();

    TEST("Add problems");
    int ov_id  = rw_add_problem(rw, "OV");
    int edit_id = rw_add_problem(rw, "Edit Distance");
    int lcs_id  = rw_add_problem(rw, "LCS");
    assert(ov_id == 0);
    assert(edit_id == 1);
    assert(lcs_id == 2);
    assert(rw->n_problems == 3);
    PASS();

    TEST("Add reductions");
    rw_add_reduction(rw, ov_id, edit_id, 1.0, 0.0, 0, "Backurs-Indyk 2015");
    rw_add_reduction(rw, ov_id, lcs_id, 1.0, 0.0, 0, "ABB 2015");
    assert(rw->n_edges == 2);
    PASS();

    TEST("Path existence");
    assert(rw_path_exists(rw, ov_id, edit_id) == 1);
    assert(rw_path_exists(rw, ov_id, lcs_id) == 1);
    assert(rw_path_exists(rw, edit_id, lcs_id) == 0);
    PASS();

    TEST("Equivalence class");
    int equiv_out[10];
    int n_equiv = rw_equivalence_class(rw, ov_id, equiv_out, 10);
    assert(n_equiv == 1); /* OV only equivalent to itself */
    assert(equiv_out[0] == ov_id);
    PASS();

    rw_free(rw);
}

/* L4: SETH theorem tests */
static void test_seth_theorems(void) {
    TEST("SETH lower bound for k=3");
    double sk3 = seth_lower_bound_for_k(3);
    assert(sk3 > 0.3 && sk3 < 0.5); /* around 0.386 */
    PASS();

    TEST("SETH lower bound for k=7");
    double sk7 = seth_lower_bound_for_k(7);
    assert(sk7 > 0.7 && sk7 < 0.85); /* around 0.775 */
    PASS();

    TEST("SETH circuit SAT bound");
    double cb = seth_circuit_sat_bound(3, 100);
    assert(cb > 0.0 && cb <= 1.0);
    PASS();

    TEST("SETH threshold for Edit Distance");
    double th = seth_threshold("Edit Distance");
    assert(th == 2.0);
    PASS();
}

/* L5: SAT algorithm tests */
static void test_sat_solvers(void) {
    TEST("SAT brute force (satisfiable)");
    SatInstance* inst = sat_create(3, 3, 3);
    int c1[] = {1, 2, 3};     /* x1 OR x2 OR x3 */
    int c2[] = {-1, -2};       /* NOT x1 OR NOT x2 */
    int c3[] = {1, -3};        /* x1 OR NOT x3 */
    sat_add_clause(inst, c1, 3);
    sat_add_clause(inst, c2, 2);
    sat_add_clause(inst, c3, 2);
    int assign[3];
    int result = sat_brute_force(inst, assign);
    assert(result == 1);
    assert(sat_check(inst, assign) == 1);
    PASS();

    TEST("SAT brute force (unsatisfiable)");
    SatInstance* inst2 = sat_create(2, 2, 2);
    int c4[] = {1};       /* x1 */
    int c5[] = {-1};      /* NOT x1 */
    sat_add_clause(inst2, c4, 1);
    sat_add_clause(inst2, c5, 1);
    int result2 = sat_brute_force(inst2, NULL);
    assert(result2 == 0);
    PASS();

    TEST("SAT DPLL (satisfiable)");
    int assign3[3];
    int result3 = sat_dpll(inst, assign3);
    assert(result3 == 1);
    assert(sat_check(inst, assign3) == 1);
    PASS();

    TEST("SAT DPLL (unsatisfiable)");
    int result4 = sat_dpll(inst2, NULL);
    assert(result4 == 0);
    PASS();

    sat_free(inst);
    sat_free(inst2);
}

/* L5: 3SUM algorithm tests */
static void test_threesum(void) {
    TEST("3SUM create/add");
    ThreeSumInstance* ts = ts3_create(10);
    assert(ts != NULL);
    ts3_add(ts, -25); ts3_add(ts, -10); ts3_add(ts, -7);
    ts3_add(ts, -3); ts3_add(ts, 2); ts3_add(ts, 4);
    ts3_add(ts, 8); ts3_add(ts, 10); ts3_add(ts, 25);
    assert(ts->n == 9);
    PASS();

    TEST("3SUM quadratic (exists)");
    ThreeSumSolution sol = ts3_quadratic(ts);
    assert(sol.found == 1); /* -25 + 0 + 25? Let's check -10 + 2 + 8 = 0 */
    int64_t s = ts->numbers[sol.i] + ts->numbers[sol.j] + ts->numbers[sol.k];
    assert(s == 0);
    PASS();

    TEST("3SUM exists check");
    assert(ts3_exists(ts) == 1);
    PASS();

    ts3_free(ts);
}

/* L5: Edit distance tests */
static void test_edit_distance(void) {
    TEST("Edit distance DP (equal strings)");
    int d1 = edit_distance_dp("hello", 5, "hello", 5);
    assert(d1 == 0);
    PASS();

    TEST("Edit distance DP (kitten/sitting)");
    int d2 = edit_distance_dp("kitten", 6, "sitting", 7);
    assert(d2 == 3); /* k->s, e->i, +g */
    PASS();

    TEST("Edit distance DP (empty string)");
    int d3 = edit_distance_dp("abc", 3, "", 0);
    assert(d3 == 3);
    PASS();

    TEST("LCS length");
    int l1 = lcs_length("ABCDGH", 6, "AEDFHR", 6);
    assert(l1 == 3); /* ADH */
    PASS();

    TEST("LCS length (no common)");
    int l2 = lcs_length("ABC", 3, "XYZ", 3);
    assert(l2 == 0);
    PASS();
}

/* L5: APSP tests */
static void test_apsp(void) {
    TEST("APSP create and Floyd-Warshall");
    ApspInstance* apsp = apsp_create(4);
    apsp_add_edge(apsp, 0, 1, 5);
    apsp_add_edge(apsp, 0, 2, 10);
    apsp_add_edge(apsp, 1, 2, 3);
    apsp_add_edge(apsp, 1, 3, 20);
    apsp_add_edge(apsp, 2, 3, 2);
    apsp_floyd_warshall(apsp);
    /* 0->1(5)->2(3)->3(2) = 10, or 0->2(10)->3(2) = 12 */
    assert(apsp->dist[0][3] == 10);
    PASS();

    TEST("APSP negative cycle check");
    assert(apsp_has_negative_cycle(apsp) == 0);
    PASS();

    apsp_free(apsp);
}

/* L5: Vector set tests */
static void test_vector_set(void) {
    TEST("Vector set create and add");
    VectorSet* vs = vs_create(10, 4);
    int v1[] = {1, 0, 0, 1};
    int v2[] = {0, 1, 1, 0};
    int v3[] = {0, 0, 0, 1};
    int id1 = vs_add_vector(vs, v1);
    int id2 = vs_add_vector(vs, v2);
    int id3 = vs_add_vector(vs, v3);
    assert(id1 == 0); assert(id2 == 1); assert(id3 == 2);
    PASS();

    TEST("Dot product");
    uint64_t dot12 = vs_dot_product(vs, 0, 1);
    assert(dot12 == 0); /* 1&0 + 0&1 + 0&1 + 1&0 = 0 (orthogonal!) */
    uint64_t dot13 = vs_dot_product(vs, 0, 2);
    assert(dot13 == 1); /* 1&0 + 0&0 + 0&0 + 1&1 = 1 */
    PASS();

    TEST("OV brute force (found)");
    OrthogonalityResult res;
    int found = ov_brute_force(vs, &res);
    assert(found == 1);
    assert(res.found == 1);
    assert(res.orthogonal == 1);
    PASS();

    vs_free(vs);
}

/* L6: Canonical problems tests */
static void test_canonical_problems(void) {
    TEST("3SUM instance count all");
    ThreeSumInstance* ts = ts3_create(6);
    ts3_add(ts, -1); ts3_add(ts, 0); ts3_add(ts, 1);
    ts3_add(ts, 2); ts3_add(ts, -2); ts3_add(ts, -3);
    int count = ts3_count_all(ts);
    assert(count >= 1);
    PASS();

    TEST("Geometric 3-collinear");
    int64_t x[] = {0, 1, 2, 3, 4};
    int64_t y[] = {0, 2, 4, 6, 1};
    int col = geom_three_collinear(x, y, 5);
    assert(col == 1); /* (0,0), (1,2), (2,4) are collinear */
    PASS();

    ts3_free(ts);
}

/* L2: Time bound tests */
static void test_time_bounds(void) {
    TEST("Time bound compare");
    FineGrainedTimeBound tb1 = fgtb_make(2.0, 0.0, 1, 0, 0.0);
    FineGrainedTimeBound tb2 = fgtb_make(3.0, 0.0, 1, 0, 0.0);
    assert(fgtb_compare(tb1, tb2) == -1);
    PASS();

    TEST("Time bound format");
    char buf[128];
    fgtb_format(tb1, buf, 128);
    assert(strstr(buf, "n") != NULL);
    PASS();
}

/* L7: Application tests */
static void test_applications(void) {
    TEST("SETH check refutation");
    int ref = seth_check_refutation(1.5, "Edit Distance", "test");
    assert(ref == 1); /* 1.5 < 2.0 would refute SETH */
    PASS();

    TEST("NSETH consistency check");
    int cons = nseth_consistency_check(0.999, 0.4);
    assert(cons == 0); /* nondet exp too small */
    PASS();
}

/* L8: Extended reduction web tests */
static void test_extended_web(void) {
    TEST("Extended web create/add");
    ExtendedReductionWeb* erw = erw_create(10, 20, 5);
    int ov = erw_add_problem(erw, "OV", CAT_STRING, 2.0, 2.0);
    int ed = erw_add_problem(erw, "Edit Distance", CAT_STRING, 2.0, 2.0);
    erw_add_edge(erw, ov, ed, 1.0, 0.0, 0, "Backurs-Indyk 2015");
    erw_set_clb(erw, ov, 2.0);
    assert(erw->n_problems == 2);
    assert(erw->n_edges == 1);
    PASS();

    TEST("SCC classes");
    int nc = erw_compute_classes(erw);
    assert(nc == 2); /* two separate SCCs */
    PASS();

    erw_free(erw);
}

int main(void) {
    printf("=== mini-conditional-lower-bounds Test Suite ===\n\n");

    printf("--- L1: Definitions ---\n");
    test_hypothesis_db();

    printf("\n--- L2: Core Concepts ---\n");
    test_time_bounds();

    printf("\n--- L3: Mathematical Structures ---\n");
    test_reduction_web();

    printf("\n--- L4: Fundamental Theorems ---\n");
    test_seth_theorems();

    printf("\n--- L5: Algorithms ---\n");
    test_sat_solvers();
    test_threesum();
    test_edit_distance();
    test_apsp();
    test_vector_set();

    printf("\n--- L6: Canonical Problems ---\n");
    test_canonical_problems();

    printf("\n--- L7: Applications ---\n");
    test_applications();

    printf("\n--- L8: Advanced Topics ---\n");
    test_extended_web();

    printf("\n=== Results: %d/%d tests passed ===\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}
