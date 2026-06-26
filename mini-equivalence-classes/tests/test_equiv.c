/* ============================================================================
 * test_equiv.c -- Comprehensive Tests for mini-equivalence-classes
 *
 * Tests cover all core APIs across L1-L8.
 * Uses assert() for self-checking tests.
 * ============================================================================ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <float.h>

#include "equiv_classes.h"
#include "subcubic.h"
#include "subquadratic.h"
#include "threesum.h"
#include "fine_grained_reduction.h"

extern void equiv_classes_init(void);

static int tests_run = 0;
static int tests_passed = 0;

#define TEST(name) do { tests_run++; printf("  TEST %2d: %-50s ", tests_run, name); } while(0)
#define PASS() do { printf("PASS\n"); tests_passed++; } while(0)
#define FAIL(msg) do { printf("FAIL: %s\n", msg); } while(0)
#define CHECK(cond) do { if (cond) PASS(); else FAIL(#cond); } while(0)

/* ---- L1: Equivalence Class Tests ---- */

static void test_equiv_class_names(void) {
    TEST("equiv_class_name valid");
    const char *n = equiv_class_name(EQUIV_CLASS_SUBCUBIC);
    CHECK(n != NULL && strlen(n) > 0);

    TEST("equiv_class_name invalid");
    n = equiv_class_name((equiv_class_id_t)999);
    CHECK(n != NULL);
}

static void test_problem_registration(void) {
    TEST("register_problem basic");
    problem_id_t pid = register_problem("TestProblem", EQUIV_CLASS_SUBCUBIC, 3.0, 100);
    CHECK(pid != (problem_id_t)-1);

    TEST("get_problem_descriptor");
    const problem_descriptor_t *desc = get_problem_descriptor(pid);
    CHECK(desc != NULL && desc->problem_id == pid);

    TEST("problem descriptor name");
    CHECK(strcmp(desc->name, "TestProblem") == 0);

    TEST("problem descriptor exponent");
    CHECK(fabs(desc->time_exponent - 3.0) < 1e-9);

    TEST("total_registered_problems");
    int32_t total = total_registered_problems();
    CHECK(total > 0);
}

static void test_equiv_class_create(void) {
    TEST("equiv_class_create valid");
    equiv_class_t *ec = equiv_class_create(EQUIV_CLASS_APSP, "TestAPSP", 3.0, 3.0);
    CHECK(ec != NULL);

    TEST("equiv_class_add_member");
    problem_id_t pid = register_problem("Member1", EQUIV_CLASS_APSP, 3.0, 10);
    bool ok = equiv_class_add_member(ec, pid);
    CHECK(ok);

    TEST("equiv_class_contains");
    CHECK(equiv_class_contains(ec, pid));

    TEST("problems_in_same_class");
    problem_id_t pid2 = register_problem("Member2", EQUIV_CLASS_APSP, 3.0, 10);
    equiv_class_add_member(ec, pid2);
    CHECK(problems_in_same_class(pid, pid2));
}

/* ---- L2: Reduction Tests ---- */

static void test_reductions(void) {
    equiv_classes_init();

    TEST("register_reduction valid");
    bool ok = register_reduction(0, 1, FGR_SUBCUBIC, 2.5, 1, "Test citation");
    CHECK(ok);

    TEST("reduction_chain_exponent self");
    double chain = reduction_chain_exponent(0, 0);
    CHECK(fabs(chain - 0.0) < 1e-9);

    TEST("is_canonical_problem");
    problem_id_t apsp_id = register_problem("APSP", EQUIV_CLASS_SUBCUBIC, 3.0, 10);
    equiv_class_t *ec = equiv_class_create(EQUIV_CLASS_SUBCUBIC, "APSP", 3.0, 3.0);
    equiv_class_add_member(ec, apsp_id);
    CHECK(is_canonical_problem(apsp_id));

    TEST("would_refute_subcubic_conjecture");
    CHECK(would_refute_subcubic_conjecture(2.5, 100));

    TEST("would_refute_ov_conjecture");
    CHECK(would_refute_ov_conjecture(1.5, 100, 100));

    TEST("would_refute_3sum_conjecture");
    CHECK(would_refute_3sum_conjecture(1.5, 100));
}

/* ---- L5: APSP / Floyd-Warshall Tests ---- */

static void test_apsp_floyd_warshall(void) {
    TEST("apsp_create");
    apsp_instance_t *apsp = apsp_create(4);
    CHECK(apsp != NULL);

    TEST("apsp_set_edge / get_edge");
    apsp_set_edge(apsp, 0, 1, 5.0);
    double w = apsp_get_edge(apsp, 0, 1);
    CHECK(fabs(w - 5.0) < 1e-9);

    TEST("apsp_get_distance self");
    double d = apsp_get_distance(apsp, 0, 0);
    CHECK(fabs(d - 0.0) < 1e-9);

    TEST("apsp_floyd_warshall triangle");
    apsp_set_edge(apsp, 1, 2, 3.0);
    apsp_set_edge(apsp, 0, 2, 10.0);
    bool ok = apsp_floyd_warshall(apsp);
    CHECK(ok);

    TEST("apsp_floyd_warshall path 0->2");
    double path = apsp_get_distance(apsp, 0, 2);
    CHECK(fabs(path - 8.0) < 1e-9);

    apsp_free(apsp);
}

/* ---- L5: Graph Metrics Tests ---- */

static void test_graph_metrics(void) {
    TEST("graph_metrics_create");
    int32_t n = 3;
    double *dist = (double *)malloc((size_t)(n * n) * sizeof(double));
    dist[0*n+0] = 0; dist[0*n+1] = 2; dist[0*n+2] = 5;
    dist[1*n+0] = 2; dist[1*n+1] = 0; dist[1*n+2] = 3;
    dist[2*n+0] = 5; dist[2*n+1] = 3; dist[2*n+2] = 0;
    graph_metrics_t *gm = graph_metrics_create(n, dist);
    CHECK(gm != NULL);

    TEST("graph_metrics_compute");
    graph_metrics_compute(gm);
    CHECK(fabs(gm->diameter - 5.0) < 1e-9);

    TEST("graph_metrics radius");
    CHECK(fabs(gm->radius - 3.0) < 1e-9);

    graph_metrics_free(gm);
    free(dist);
}

/* ---- L5: Min-Plus Product Tests ---- */

static void test_min_plus_product(void) {
    TEST("min_plus_create");
    min_plus_product_t *mp = min_plus_create(2);
    CHECK(mp != NULL);

    TEST("min_plus_naive");
    mp->A[0*2+0] = 0; mp->A[0*2+1] = 3;
    mp->A[1*2+0] = 5; mp->A[1*2+1] = 0;
    mp->B[0*2+0] = 0; mp->B[0*2+1] = 1;
    mp->B[1*2+0] = 2; mp->B[1*2+1] = 0;
    min_plus_naive(mp);
    CHECK(fabs(mp->C[0*2+0] - 0.0) < 1e-9);
    CHECK(fabs(mp->C[0*2+1] - 1.0) < 1e-9);

    TEST("min_plus_verify");
    min_plus_product_t *mp2 = min_plus_create(2);
    memcpy(mp2->A, mp->A, 4*sizeof(double));
    memcpy(mp2->B, mp->B, 4*sizeof(double));
    min_plus_naive(mp2);
    CHECK(min_plus_verify(mp, mp2, 1e-6));

    min_plus_free(mp);
    min_plus_free(mp2);
}

/* ---- L5: Orthogonal Vectors Tests ---- */

static void test_ov(void) {
    TEST("ov_create");
    ov_instance_t *ov = ov_create(5, 8);
    CHECK(ov != NULL);

    TEST("ov_set/get_bit");
    ov_set_bit(ov, 0, 3, true, true);
    CHECK(ov_get_bit(ov, 0, 3, true) == true);

    TEST("ov_brute_force no pair");
    for (int32_t k = 0; k < 5; k++) ov_set_bit(ov, k, 0, true, true);
    int32_t ai, bi;
    bool found = ov_brute_force(ov, &ai, &bi);
    CHECK(found == false);

    ov_free(ov);
}

/* ---- L5: Edit Distance Tests ---- */

static void test_edit_distance(void) {
    TEST("edit_distance_levenshtein kitten/sitting");
    int32_t d = edit_distance_levenshtein("kitten", 6, "sitting", 7);
    CHECK(d == 3);

    TEST("edit_distance_levenshtein same");
    d = edit_distance_levenshtein("abc", 3, "abc", 3);
    CHECK(d == 0);

    TEST("edit_distance_levenshtein empty");
    d = edit_distance_levenshtein("abc", 3, "", 0);
    CHECK(d == 3);

    TEST("edit_distance_dp full");
    edit_distance_instance_t *ed = edit_distance_create(3, 3);
    strncpy(ed->x, "abc", 3);
    strncpy(ed->y, "abd", 3);
    d = edit_distance_dp(ed);
    CHECK(d == 1);
    edit_distance_free(ed);
}

/* ---- L5: LCS Tests ---- */

static void test_lcs(void) {
    TEST("lcs_dp ABCD/ACBD");
    lcs_instance_t *lcs = lcs_create(4, 4);
    strncpy(lcs->x, "ABCD", 4);
    strncpy(lcs->y, "ACBD", 4);
    int32_t len = lcs_dp(lcs);
    CHECK(len == 3);

    TEST("lcs_hirschberg");
    lcs_instance_t *lcs2 = lcs_create(4, 4);
    strncpy(lcs2->x, "ABCD", 4);
    strncpy(lcs2->y, "ACBD", 4);
    len = lcs_hirschberg(lcs2);
    CHECK(len == 3);
    lcs_free(lcs);
    lcs_free(lcs2);
}

/* ---- L5: DTW Tests ---- */

static void test_dtw(void) {
    TEST("dtw_create");
    dtw_instance_t *dtw = dtw_create(3, 3);
    CHECK(dtw != NULL);

    TEST("dtw_standard");
    dtw->x[0] = 1.0; dtw->x[1] = 2.0; dtw->x[2] = 3.0;
    dtw->y[0] = 1.0; dtw->y[1] = 2.0; dtw->y[2] = 3.0;
    double d = dtw_standard(dtw);
    CHECK(fabs(d - 0.0) < 1e-9);
    dtw_free(dtw);
}

/* ---- L5: 3SUM Tests ---- */

static void test_3sum(void) {
    TEST("threesum_create");
    threesum_instance_t *ts = threesum_create(5);
    CHECK(ts != NULL);

    TEST("threesum_quadratic found");
    ts->values[0] = -10; ts->values[1] = 3;
    ts->values[2] = 5; ts->values[3] = 7; ts->values[4] = 2;
    bool found = threesum_quadratic(ts);
    CHECK(found == true);

    TEST("threesum_quadratic not found");
    threesum_instance_t *ts2 = threesum_create(3);
    ts2->values[0] = 1; ts2->values[1] = 2; ts2->values[2] = 3;
    found = threesum_quadratic(ts2);
    CHECK(found == false);

    threesum_free(ts);
    threesum_free(ts2);
}

/* ---- L7: Application Tests ---- */

static void test_applications(void) {
    TEST("dna_sequence_alignment");
    int32_t score = dna_sequence_alignment("ACGT", 4, "AGCT", 4, 2, -1, -1);
    CHECK(score >= 0);

    TEST("mars_rover_collinearity check");
    double ox[] = {0.0, 1.0, 2.0, 3.0, 4.0};
    double oy[] = {0.0, 1.0, 2.0, 3.0, 4.0};
    bool collinear = mars_rover_collinearity_obstacle_check(ox, oy, 5);
    CHECK(collinear == true);
}

/* ---- L8: Fine-Grained Reduction Graph Tests ---- */

static void test_fg_graph(void) {
    TEST("fg_graph_create");
    fg_reduction_graph_t *g = fg_graph_create(10, 20);
    CHECK(g != NULL);

    TEST("fg_graph_add_node");
    fg_node_t *n = fg_graph_add_node(g, "Node1", EQUIV_CLASS_SUBCUBIC, 3.0, 3.0);
    CHECK(n != NULL);

    TEST("fg_graph_find_node");
    fg_node_t *f = fg_graph_find_node(g, "Node1");
    CHECK(f != NULL && f == n);

    TEST("fg_graph_add_edge");
    fg_node_t *n2 = fg_graph_add_node(g, "Node2", EQUIV_CLASS_SUBCUBIC, 2.5, 3.0);
    fg_edge_t *e = fg_graph_add_edge(g, n->id, n2->id, FGR_SUBCUBIC, 0.5, 1, true, "test");
    CHECK(e != NULL);

    TEST("fg_shortest_reduction_chain");
    problem_id_t path[10];
    int32_t plen;
    double chain = fg_shortest_reduction_chain(g, 0, 1, path, &plen, 10);
    CHECK(chain >= 0);

    TEST("fg_transitive_closure");
    double **tc = fg_transitive_closure(g);
    CHECK(tc != NULL);
    fg_transitive_closure_free(tc, 2);

    TEST("fg_find_equivalence_classes");
    int32_t nclasses;
    problem_id_t *rep = fg_find_equivalence_classes(g, &nclasses);
    CHECK(rep != NULL && nclasses > 0);
    free(rep);

    fg_graph_free(g);
}

/* ---- Main Test Runner ---- */

int main(void) {
    printf("=== mini-equivalence-classes: Unit Tests ===\n\n");

    printf("-- L1: Definitions & Equivalence Classes --\n");
    test_equiv_class_names();
    test_problem_registration();
    test_equiv_class_create();

    printf("\n-- L2: Reductions & Conjectures --\n");
    test_reductions();

    printf("\n-- L5: APSP & Graph Metrics --\n");
    test_apsp_floyd_warshall();
    test_graph_metrics();
    test_min_plus_product();

    printf("\n-- L5: Subquadratic Algorithms --\n");
    test_ov();
    test_edit_distance();
    test_lcs();
    test_dtw();

    printf("\n-- L5: 3SUM Algorithms --\n");
    test_3sum();

    printf("\n-- L7: Applications --\n");
    test_applications();

    printf("\n-- L8: Fine-Grained Reduction Graph --\n");
    test_fg_graph();

    printf("\n=== Results: %d/%d tests passed ===\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}