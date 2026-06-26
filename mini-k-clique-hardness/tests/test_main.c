/**
 * test_main.c — Tests for mini-k-clique-hardness module
 *
 * Tests graph operations, clique verification, brute-force search,
 * ETH/SETH lower bound derivations, and parameterized complexity.
 *
 * All tests use standard assert() with descriptive messages.
 *
 * Knowledge: L1-L6 coverage verification via assert-based testing
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>

#include "../include/kclique_types.h"
#include "../include/kclique_core.h"
#include "../include/kclique_eth.h"
#include "../include/kclique_algorithm.h"
#include "../include/kclique_reduction.h"
#include "../include/kclique_graph.h"

/* Test counter */
static int tests_run = 0;
static int tests_passed = 0;
static int tests_failed = 0;

#define TEST(name) do { \
    tests_run++; \
    printf("  TEST %s ... ", name); \
} while(0)

#define PASS() do { \
    tests_passed++; \
    printf("PASS\n"); \
} while(0)

#define FAIL(msg) do { \
    tests_failed++; \
    printf("FAIL: %s\n", msg); \
} while(0)

#define CHECK(cond, msg) do { \
    if (!(cond)) { FAIL(msg); } else { PASS(); } \
} while(0)

/* ================================================================
 * L1: Graph Creation and Basic Operations
 * ================================================================ */

static void test_graph_create(void) {
    TEST("graph_create small");
    graph_t *g = graph_create(5);
    assert(g != NULL);
    assert(g->n == 5);
    assert(g->m == 0);
    assert(!g->directed);
    assert(g->simple);
    graph_destroy(g);
    PASS();

    TEST("graph_create zero (invalid)");
    g = graph_create(0);
    assert(g == NULL);
    PASS();

    TEST("graph_create max");
    g = graph_create(100);
    assert(g != NULL);
    assert(g->n == 100);
    graph_destroy(g);
    PASS();
}

static void test_graph_add_edge(void) {
    TEST("graph_add_edge basic");
    graph_t *g = graph_create(5);
    assert(g != NULL);
    bool r = graph_add_edge(g, 0, 1);
    assert(r == true);
    assert(g->m == 1);
    assert(graph_has_edge(g, 0, 1));
    assert(graph_has_edge(g, 1, 0)); /* undirected */
    PASS();

    TEST("graph_add_edge duplicate");
    r = graph_add_edge(g, 0, 1);
    assert(r == false); /* already exists */
    assert(g->m == 1);
    PASS();

    TEST("graph_add_edge self-loop (rejected)");
    r = graph_add_edge(g, 0, 0);
    assert(r == false);
    PASS();

    TEST("graph_add_edge out of bounds");
    r = graph_add_edge(g, -1, 0);
    assert(r == false);
    r = graph_add_edge(g, 0, 5);
    assert(r == false);
    PASS();

    graph_destroy(g);
}

static void test_graph_degree_neighbors(void) {
    TEST("graph_degree");
    graph_t *g = graph_create(5);
    graph_add_edge(g, 0, 1);
    graph_add_edge(g, 0, 2);
    graph_add_edge(g, 0, 3);
    assert(graph_degree(g, 0) == 3);
    assert(graph_degree(g, 1) == 1);
    assert(graph_degree(g, 4) == 0);
    PASS();

    TEST("graph_neighbors");
    int32_t count;
    int32_t *nbrs = graph_neighbors(g, 0, &count);
    assert(nbrs != NULL);
    assert(count == 3);
    free(nbrs);

    nbrs = graph_neighbors(g, 4, &count);
    assert(nbrs == NULL);
    assert(count == 0);
    PASS();

    graph_destroy(g);
}

static void test_graph_complement(void) {
    TEST("graph_complement");
    graph_t *g = graph_create(4);
    graph_add_edge(g, 0, 1);
    graph_add_edge(g, 0, 2);
    /* Edges: (0,1), (0,2), (1,0), (2,0) */
    /* Complement should have: (0,3), (1,2), (1,3), (2,3) */

    graph_t *comp = graph_complement(g);
    assert(comp != NULL);
    assert(comp->n == 4);
    /* Original edges should NOT be in complement */
    assert(!graph_has_edge(comp, 0, 1));
    assert(!graph_has_edge(comp, 0, 2));
    /* Non-edges should become edges in complement */
    assert(graph_has_edge(comp, 0, 3));
    assert(graph_has_edge(comp, 1, 2));
    assert(graph_has_edge(comp, 1, 3));
    assert(graph_has_edge(comp, 2, 3));
    /* Diagonal never has edges */
    assert(!graph_has_edge(comp, 0, 0));

    graph_destroy(g);
    graph_destroy(comp);
    PASS();
}

/* ================================================================
 * L1/L5: k-Clique Verification and Search
 * ================================================================ */

static void test_is_clique(void) {
    TEST("is_clique triangle");
    graph_t *g = graph_create(4);
    graph_add_edge(g, 0, 1);
    graph_add_edge(g, 0, 2);
    graph_add_edge(g, 1, 2);
    /* Triangle on {0,1,2} */

    int32_t tri[3] = {0, 1, 2};
    assert(is_clique(g, tri, 3));
    assert(is_k_clique(g, tri, 3));

    int32_t not_clique[3] = {0, 1, 3}; /* no edge (2,3) or (1,3) */
    assert(!is_clique(g, not_clique, 3));

    int32_t single[1] = {0};
    assert(is_clique(g, single, 1)); /* singleton always clique */
    PASS();

    graph_destroy(g);
}

static void test_find_k_clique_brute(void) {
    TEST("find_k_clique_brute triangle");
    graph_t *g = graph_create(5);
    graph_add_edge(g, 0, 1);
    graph_add_edge(g, 0, 2);
    graph_add_edge(g, 1, 2);
    graph_add_edge(g, 3, 4);

    clique_t result = {0};
    bool found = find_k_clique_brute(g, 3, &result);
    assert(found);
    assert(result.found);
    assert(result.k == 3);
    assert(result.vertices != NULL);
    /* Verify it's actually a clique */
    assert(is_k_clique(g, result.vertices, 3));
    clique_free(&result);
    PASS();

    TEST("find_k_clique_brute none");
    found = find_k_clique_brute(g, 4, &result);
    assert(!found);
    PASS();

    TEST("find_k_clique_brute k=1");
    found = find_k_clique_brute(g, 1, &result);
    assert(found);
    assert(result.k == 1);
    clique_free(&result);
    PASS();

    graph_destroy(g);
}

static void test_count_k_cliques(void) {
    TEST("count_k_cliques triangle");
    graph_t *g = graph_create(4);
    graph_add_edge(g, 0, 1);
    graph_add_edge(g, 0, 2);
    graph_add_edge(g, 1, 2);
    g->m = 3;

    int64_t cnt = count_k_cliques(g, 3);
    assert(cnt == 1); /* exactly one triangle */

    cnt = count_k_cliques(g, 2);
    assert(cnt == 3); /* 3 edges = 3 2-cliques */

    cnt = count_k_cliques(g, 1);
    assert(cnt == 4); /* 4 vertices = 4 1-cliques */
    PASS();

    graph_destroy(g);
}

/* ================================================================
 * L2/L4: ETH/SETH Lower Bounds
 * ================================================================ */

static void test_eth_lower_bound(void) {
    TEST("eth_kclique_lower_bound");
    double alpha = eth_kclique_lower_bound(100, 5, 1.0);
    /* Under ETH, k-Clique with k=5 on n=100 should have alpha > 0 */
    assert(alpha > 0.0);
    assert(alpha <= 5.0); /* cannot exceed k */
    PASS();

    TEST("eth_s_k_for_clause_width");
    double s_k = eth_s_k_for_clause_width(3);
    assert(s_k > 0.0);
    s_k = eth_s_k_for_clause_width(10);
    assert(s_k > 0.0);
    PASS();
}

static void test_seth_lower_bound(void) {
    TEST("seth_kclique_lower_bound");
    double lower = seth_kclique_lower_bound(1000, 4, 0.01);
    assert(lower >= 4.0 - 0.01); /* k - eps */
    assert(lower <= 4.0);
    PASS();

    TEST("seth_tightness_gap");
    double gap = seth_tightness_gap(1000000, 5);
    assert(gap >= 0.0 && gap < 1.0);
    PASS();
}

static void test_conditional_lower_bound(void) {
    TEST("conditional_lower_bound create/free");
    conditional_lower_bound_t clb = make_conditional_lower_bound(
        "ETH", "k-Clique", 5.0, false, false);
    assert(clb.conjecture != NULL);
    assert(clb.problem != NULL);
    assert(clb.time_exponent == 5.0);
    assert(!clb.exponential);
    assert(!clb.tight);
    free_clb(&clb);
    PASS();

    TEST("compare_clb_hardness");
    conditional_lower_bound_t clb1 = make_conditional_lower_bound(
        "ETH", "k-Clique", 3.0, false, false);
    conditional_lower_bound_t clb2 = make_conditional_lower_bound(
        "SETH", "k-Clique", 4.0, false, true);
    assert(compare_clb_hardness(&clb1, &clb2) < 0); /* clb2 is stronger */
    free_clb(&clb1); free_clb(&clb2);

    conditional_lower_bound_t exp_clb = make_conditional_lower_bound(
        "EXP", "Clique", 0.0, true, false);
    assert(compare_clb_hardness(&exp_clb, &clb1) > 0); /* exponential stronger */
    free_clb(&exp_clb);
    PASS();
}

/* ================================================================
 * L5: Bron-Kerbosch and Color-Coding
 * ================================================================ */

static void test_bron_kerbosch(void) {
    TEST("max_clique_bron_kerbosch");
    graph_t *g = graph_create(6);
    /* Triangle {0,1,2} + edge {3,4} + isolated {5} */
    graph_add_edge(g, 0, 1); graph_add_edge(g, 0, 2);
    graph_add_edge(g, 1, 2);
    graph_add_edge(g, 3, 4);
    g->m = 4;

    clique_t result = {0};
    int32_t max_size = max_clique_bron_kerbosch(g, &result);
    assert(max_size == 3);
    assert(result.found);
    assert(result.k == 3);
    assert(is_k_clique(g, result.vertices, 3));
    clique_free(&result);
    PASS();

    graph_destroy(g);
}

static void test_color_coding(void) {
    TEST("color_coding_solve triangle");
    graph_t *g = graph_create(6);
    graph_add_edge(g, 0, 1); graph_add_edge(g, 0, 2);
    graph_add_edge(g, 1, 2); /* triangle on {0,1,2} */
    graph_add_edge(g, 3, 4); /* just an edge */
    g->m = 4;

    clique_t result = {0};
    bool found = color_coding_solve(g, 3, 200, 42, &result);
    assert(found);
    assert(result.found);
    assert(result.k == 3);
    assert(is_k_clique(g, result.vertices, 3));
    clique_free(&result);
    PASS();

    TEST("color_coding_solve no 4-clique");
    found = color_coding_solve(g, 4, 50, 42, &result);
    assert(!found);
    PASS();

    graph_destroy(g);
}

static void test_matrix_triangle(void) {
    TEST("matrix_triangle_detect");
    graph_t *g = graph_create(5);
    graph_add_edge(g, 0, 1); graph_add_edge(g, 0, 2);
    graph_add_edge(g, 1, 2); /* triangle {0,1,2} */
    graph_add_edge(g, 2, 3); /* edge, not triangle */
    graph_add_edge(g, 3, 4);
    g->m = 5;

    clique_t result = {0};
    bool found = matrix_triangle_detect(g, &result);
    assert(found);
    assert(result.k == 3);
    assert(is_k_clique(g, result.vertices, 3));
    clique_free(&result);
    PASS();

    TEST("matrix_4clique_detect none");
    found = matrix_4clique_detect(g, &result);
    assert(!found);
    PASS();

    graph_destroy(g);
}

/* ================================================================
 * L4/L6: Reductions
 * ================================================================ */

static void test_reduction_clique_is(void) {
    TEST("kclique_to_kindependentset");
    /* G: triangle on {0,1,2}. Complement: empty triangle -> 3 isolated vertices.
       A 3-clique in G exists (true). In complement: a 3-independent set exists.
       Complement of G has no edges among {0,1,2}, so {0,1,2} is an independent set
       of size 3 in complement(G). */
    graph_t *g = graph_create(4);
    graph_add_edge(g, 0, 1);
    graph_add_edge(g, 0, 2);
    graph_add_edge(g, 1, 2);
    g->m = 3;

    graph_t *comp = kclique_to_kindependentset(g);
    assert(comp != NULL);
    /* In complement, {0,1,2} should be an independent set (no edges among them) */
    assert(!graph_has_edge(comp, 0, 1));
    assert(!graph_has_edge(comp, 0, 2));
    assert(!graph_has_edge(comp, 1, 2));
    graph_destroy(g);
    graph_destroy(comp);
    PASS();
}

static void test_reduction_clique_vc(void) {
    TEST("kclique_to_vertexcover");
    graph_t *g = graph_create(4);
    graph_add_edge(g, 0, 1);
    graph_add_edge(g, 0, 2);
    graph_add_edge(g, 1, 2); /* 3-clique on {0,1,2} */
    graph_add_edge(g, 3, 0); /* extra edge */
    g->m = 4;

    int32_t vc_k;
    graph_t *vc_graph = kclique_to_vertexcover(g, 3, &vc_k);
    assert(vc_graph != NULL);
    assert(vc_k == 1); /* n - k = 4 - 3 = 1 */
    graph_destroy(g);
    graph_destroy(vc_graph);
    PASS();
}

/* ================================================================
 * L3/L5: Graph Properties
 * ================================================================ */

static void test_graph_properties(void) {
    TEST("graph_degeneracy");
    graph_t *g = graph_create(5);
    graph_add_edge(g, 0, 1); graph_add_edge(g, 0, 2);
    graph_add_edge(g, 1, 2); /* triangle */
    graph_add_edge(g, 3, 4); /* edge */
    g->m = 4;

    int32_t *ordering;
    int32_t d = graph_degeneracy(g, &ordering);
    assert(d >= 1); /* at least degeneracy 1 for any non-empty graph */
    free(ordering);
    PASS();

    graph_destroy(g);

    TEST("graph_is_triangle_free");
    graph_t *g2 = graph_create(5);
    graph_add_edge(g2, 0, 1);
    graph_add_edge(g2, 1, 2);
    graph_add_edge(g2, 2, 3);
    graph_add_edge(g2, 3, 4);
    /* Path graph: no triangles */
    g2->m = 4;
    assert(graph_is_triangle_free(g2));
    graph_destroy(g2);
    PASS();
}

static void test_turan_bound(void) {
    TEST("turan_bound_edges");
    /* Turan bound for K_4-free graph on 10 vertices:
       ex(10, 4) = (1 - 1/3) * 10*9/2 = 2/3 * 45 = 30 */
    int64_t bound = turan_bound_edges(10, 3); /* r = k-1 = 3 for K4-free */
    assert(bound > 0);
    PASS();

    TEST("turan_guarantees_clique");
    graph_t *g = graph_generate_complete(5);
    assert(turan_guarantees_clique(g, 5)); /* K5 guaranteed to have K5 */
    graph_destroy(g);
    PASS();
}

/* ================================================================
 * Graph Generation
 * ================================================================ */

static void test_graph_generation(void) {
    TEST("graph_generate_erdos_renyi");
    er_graph_params_t params = {50, 0.3, 42, false, false, 0};
    graph_t *g = graph_generate_erdos_renyi(&params);
    assert(g != NULL);
    assert(g->n == 50);
    assert(g->m > 0); /* with p=0.3, should have some edges */
    graph_destroy(g);
    PASS();

    TEST("graph_generate_complete");
    graph_t *k5 = graph_generate_complete(5);
    assert(k5 != NULL);
    assert(k5->n == 5);
    /* K5 has 5*4/2 = 10 edges (but m tracks added calls, each call to graph_add_edge
       from graph_generate_complete increments m via add_edge which returns true) */
    /* Let's just check that all pairs have edges */
    for (int32_t i = 0; i < 5; i++) {
        for (int32_t j = i + 1; j < 5; j++) {
            assert(graph_has_edge(k5, i, j));
        }
    }
    graph_destroy(k5);
    PASS();

    TEST("graph_density");
    graph_t *g3 = graph_generate_empty(10);
    assert(graph_density(g3) == 0.0);
    graph_t *g4 = graph_generate_complete(10);
    assert(graph_density(g4) == 1.0);
    graph_destroy(g3);
    graph_destroy(g4);
    PASS();
}

/* ================================================================
 * Parameterized Complexity
 * ================================================================ */

static void test_parameterized(void) {
    TEST("fpt_param_blowup_kclique");
    int32_t blowup = fpt_param_blowup_kclique(3);
    assert(blowup > 0);
    PASS();

    TEST("is_polynomial_kernel false for large k");
    /* k-Clique is NOT expected to have polynomial kernel */
    /* is_polynomial_kernel checks size <= k^10 */
    assert(!is_polynomial_kernel(10000000, 3));
    PASS();
}

/* ================================================================
 * Utility: Common Neighborhood
 * ================================================================ */

static void test_common_neighborhood(void) {
    TEST("common_neighborhood");
    graph_t *g = graph_create(6);
    /* Star from vertex 0 to 1,2,3. Edge (1,2). */
    graph_add_edge(g, 0, 1);
    graph_add_edge(g, 0, 2);
    graph_add_edge(g, 0, 3);
    graph_add_edge(g, 1, 2);
    graph_add_edge(g, 4, 5);
    g->m = 5;

    int32_t set[2] = {0, 1};
    int32_t count;
    int32_t *cn = common_neighborhood(g, set, 2, &count);
    /* Common neighbors of {0,1}: vertices adjacent to both.
       N(0) = {1,2,3}, N(1) = {0,2}, common = {2} */
    assert(cn != NULL);
    assert(count == 1);
    assert(cn[0] == 2);
    free(cn);
    PASS();

    graph_destroy(g);
}

/* ================================================================
 * Main Test Runner
 * ================================================================ */

int main(void) {
    printf("=== mini-k-clique-hardness Test Suite ===\n\n");

    printf("--- L1: Graph Operations ---\n");
    test_graph_create();
    test_graph_add_edge();
    test_graph_degree_neighbors();
    test_graph_complement();

    printf("\n--- L1/L5: k-Clique Verification and Search ---\n");
    test_is_clique();
    test_find_k_clique_brute();
    test_count_k_cliques();

    printf("\n--- L2/L4: ETH/SETH Lower Bounds ---\n");
    test_eth_lower_bound();
    test_seth_lower_bound();
    test_conditional_lower_bound();

    printf("\n--- L5: Algorithms ---\n");
    test_bron_kerbosch();
    test_color_coding();
    test_matrix_triangle();

    printf("\n--- L4/L6: Reductions ---\n");
    test_reduction_clique_is();
    test_reduction_clique_vc();

    printf("\n--- L3/L5: Graph Properties ---\n");
    test_graph_properties();
    test_turan_bound();
    test_common_neighborhood();

    printf("\n--- Graph Generation ---\n");
    test_graph_generation();

    printf("\n--- Parameterized Complexity ---\n");
    test_parameterized();

    printf("\n========================================\n");
    printf("Tests run:    %d\n", tests_run);
    printf("Tests passed: %d\n", tests_passed);
    printf("Tests failed: %d\n", tests_failed);

    if (tests_failed > 0) {
        printf("\n*** SOME TESTS FAILED ***\n");
        return 1;
    }

    printf("\n*** ALL TESTS PASSED ***\n");
    return 0;
}
