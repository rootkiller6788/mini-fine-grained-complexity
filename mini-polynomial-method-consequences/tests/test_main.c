#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "poly_method.h"

/* ============================================================================
 * test_main.c - Comprehensive Tests for Polynomial Method Module
 *
 * Tests cover: L1 definitions, L3 mathematical structures,
 * L4 fundamental laws, L5 algorithms, L6 canonical problems.
 * ============================================================================ */

static int tests_passed = 0;
static int tests_total = 0;
#define TEST(name) do { tests_total++; printf("  TEST %d: %s ... ", tests_total, name); } while(0)
#define PASS() do { printf("PASS\n"); tests_passed++; } while(0)
#define FAIL(msg) do { printf("FAIL: %s\n", msg); } while(0)

/* ---- L1: Type Definitions and Creation/Destruction ---- */
static void test_l1_definitions(void) {
    TEST("poly_create/destroy basic");
    polynomial_t *p = poly_create(4, 16, false);
    assert(p != NULL);
    assert(p->num_terms == 0);
    assert(p->num_vars == 4);
    assert(p->degree == 0);
    assert(!p->over_gf2);
    poly_destroy(p);
    PASS();

    TEST("poly_create GF(2) mode");
    polynomial_t *p2 = poly_create(8, 32, true);
    assert(p2 != NULL);
    assert(p2->over_gf2);
    poly_destroy(p2);
    PASS();

    TEST("poly_clone");
    polynomial_t *orig = poly_create(3, 10, false);
    poly_add_term(orig, 1.0, 0x1);
    poly_add_term(orig, 2.0, 0x3);
    polynomial_t *clone = poly_clone(orig);
    assert(clone != NULL);
    assert(clone->num_terms == orig->num_terms);
    assert(fabs(poly_evaluate_bool(clone, 0x1) - poly_evaluate_bool(orig, 0x1)) < 1e-9);
    poly_destroy(orig);
    poly_destroy(clone);
    PASS();

    TEST("ov_create/destroy");
    ov_instance_t *ov = ov_create(10, 8);
    assert(ov != NULL);
    assert(ov->num_vectors == 10);
    assert(ov->dimension == 8);
    ov_destroy(ov);
    PASS();

    TEST("williams_cnf_create/destroy");
    williams_cnf_t *f = williams_cnf_create(5, 3, 3);
    assert(f != NULL);
    assert(f->num_vars == 5);
    assert(f->clause_width == 3);
    williams_cnf_destroy(f);
    PASS();
}

/* ---- L3: GF(p) Arithmetic ---- */
static void test_gf_arithmetic(void) {
    TEST("gf_p_add basic");
    assert(gf_p_add(3, 5, 7) == 1);
    assert(gf_p_add(6, 4, 7) == 3);
    PASS();

    TEST("gf_p_mul basic");
    assert(gf_p_mul(3, 5, 7) == 1);
    assert(gf_p_mul(2, 4, 7) == 1);
    assert(gf_p_mul(0, 5, 7) == 0);
    PASS();

    TEST("gf_p_pow Fermat's little theorem");
    assert(gf_p_pow(3, 6, 7) == 1); /* 3^6 ? 1 mod 7 */
    assert(gf_p_pow(2, 10, 11) == 1); /* 2^10 ? 1 mod 11 */
    PASS();

    TEST("gf_p_neg");
    assert(gf_p_neg(3, 7) == 4);
    assert(gf_p_neg(0, 7) == 0);
    PASS();

    TEST("gf_p_inv");
    assert(gf_p_inv(3, 7) == 5);  /* 3*5=15?1 mod 7 */
    assert(gf_p_inv(1, 7) == 1);
    assert(gf_p_inv(0, 7) == 0);  /* no inverse */
    assert(gf_p_inv(7, 11) == 8); /* 7*8=56?1 mod 11 */
    PASS();
}

/* ---- L3: Polynomial Operations ---- */
static void test_polynomial_operations(void) {
    TEST("poly_add_term basic");
    polynomial_t *p = poly_create(3, 10, false);
    poly_add_term(p, 3.0, 0x1);   /* 3*x0 */
    poly_add_term(p, 5.0, 0x2);   /* 5*x1 */
    assert(p->num_terms == 2);
    assert(p->degree == 1);
    poly_destroy(p);
    PASS();

    TEST("poly_add_term combining like terms");
    polynomial_t *p2 = poly_create(2, 10, false);
    poly_add_term(p2, 3.0, 0x1);
    poly_add_term(p2, 4.0, 0x1);  /* should combine to 7*x0 */
    assert(p2->num_terms == 1);
    assert(fabs(p2->terms[0].coeff - 7.0) < 1e-9);
    poly_destroy(p2);
    PASS();

    TEST("poly_add_term GF(2) XOR cancellation");
    polynomial_t *p3 = poly_create(2, 10, true);
    poly_add_term(p3, 1.0, 0x1);
    poly_add_term(p3, 1.0, 0x1);  /* XOR: 1 XOR 1 = 0, term removed */
    assert(p3->num_terms == 0);
    poly_destroy(p3);
    PASS();

    TEST("poly_evaluate_bool simple");
    polynomial_t *p4 = poly_create(2, 10, false);
    poly_add_term(p4, 1.0, 0x0);  /* constant 1 */
    poly_add_term(p4, 2.0, 0x1);  /* 2*x0 */
    poly_add_term(p4, 3.0, 0x2);  /* 3*x1 */
    /* p4 = 1 + 2*x0 + 3*x1 */
    assert(fabs(poly_evaluate_bool(p4, 0x0) - 1.0) < 1e-9);
    assert(fabs(poly_evaluate_bool(p4, 0x1) - 3.0) < 1e-9);  /* 1+2 */
    assert(fabs(poly_evaluate_bool(p4, 0x2) - 4.0) < 1e-9);  /* 1+3 */
    assert(fabs(poly_evaluate_bool(p4, 0x3) - 6.0) < 1e-9);  /* 1+2+3 */
    poly_destroy(p4);
    PASS();

    TEST("poly_evaluate_gf2");
    polynomial_t *p5 = poly_create(2, 10, true);
    poly_add_term(p5, 1.0, 0x1);  /* x0 */
    poly_add_term(p5, 1.0, 0x2);  /* x1 */
    /* p5 = x0 XOR x1 */
    assert(poly_evaluate_gf2(p5, 0x0) == 0);
    assert(poly_evaluate_gf2(p5, 0x1) == 1);
    assert(poly_evaluate_gf2(p5, 0x2) == 1);
    assert(poly_evaluate_gf2(p5, 0x3) == 0);  /* 1 XOR 1 = 0 */
    poly_destroy(p5);
    PASS();

    TEST("poly_add");
    polynomial_t *a = poly_create(3, 10, false);
    polynomial_t *b = poly_create(3, 10, false);
    poly_add_term(a, 1.0, 0x1);
    poly_add_term(b, 2.0, 0x2);
    polynomial_t *c = poly_add(a, b);
    assert(c->num_terms == 2);
    assert(fabs(poly_evaluate_bool(c, 0x1) - 1.0) < 1e-9);
    assert(fabs(poly_evaluate_bool(c, 0x2) - 2.0) < 1e-9);
    poly_destroy(a); poly_destroy(b); poly_destroy(c);
    PASS();

    TEST("poly_mul");
    polynomial_t *m1 = poly_create(2, 10, false);
    polynomial_t *m2 = poly_create(2, 10, false);
    poly_add_term(m1, 1.0, 0x1);  /* x0 */
    poly_add_term(m2, 1.0, 0x2);  /* x1 */
    polynomial_t *prod = poly_mul(m1, m2);  /* x0 * x1 */
    assert(prod->num_terms >= 1);  /* may have combined terms */
    assert(fabs(poly_evaluate_bool(prod, 0x3) - 1.0) < 1e-9);
    assert(fabs(poly_evaluate_bool(prod, 0x1) - 0.0) < 1e-9);
    poly_destroy(m1); poly_destroy(m2); poly_destroy(prod);
    PASS();

    TEST("poly_make_multilinear");
    polynomial_t *ml = poly_create(2, 10, true);
    poly_add_term(ml, 1.0, 0x1);
    poly_add_term(ml, 1.0, 0x3);  /* x0*x1 */
    poly_make_multilinear(ml);
    assert(poly_is_multilinear(ml));
    poly_destroy(ml);
    PASS();
}

/* ---- L4: Polynomial Representations of Boolean Functions ---- */
static void test_polynomial_representations(void) {
    TEST("truth_table_to_poly_gf2 AND function");
    uint8_t tt_and[8] = {0,0,0,0,0,0,0,1}; /* 3-var AND */
    polynomial_t *p = truth_table_to_poly_gf2(tt_and, 3);
    assert(p != NULL);
    /* AND = x0*x1*x2, should be a single monomial of degree 3 */
    assert(poly_degree(p) == 3);
    /* Verify evaluation at all points */
    for (int32_t x = 0; x < 8; x++) {
        assert(poly_evaluate_gf2(p, (uint64_t)x) == tt_and[x]);
    }
    poly_destroy(p);
    PASS();

    TEST("truth_table_to_poly_gf2 PARITY function");
    uint8_t tt_par[8] = {0,1,1,0,1,0,0,1}; /* 3-var XOR */
    polynomial_t *p2 = truth_table_to_poly_gf2(tt_par, 3);
    assert(p2 != NULL);
    /* PARITY = x0 + x1 + x2 (mod 2), degree 1 */
    assert(poly_degree(p2) == 1);
    for (int32_t x = 0; x < 8; x++) {
        assert(poly_evaluate_gf2(p2, (uint64_t)x) == tt_par[x]);
    }
    poly_destroy(p2);
    PASS();

    TEST("truth_table_to_poly_fourier AND function");
    uint8_t tt_and2[8] = {0,0,0,0,0,0,0,1};
    polynomial_t *pf = truth_table_to_poly_fourier(tt_and2, 3);
    assert(pf != NULL);
    assert(pf->num_terms > 0);
    poly_destroy(pf);
    PASS();

    TEST("fourier_coefficient constant function");
    uint8_t tt_const[4] = {1,1,1,1}; /* n=2, constant 1 */
    /* f_hat(empty) should be 1 */
    double f0 = fourier_coefficient(tt_const, 2, 0ULL);
    assert(fabs(f0 - 1.0) < 1e-9);
    /* Other coefficients should be 0 */
    double f1 = fourier_coefficient(tt_const, 2, 1ULL);
    assert(fabs(f1) < 1e-9);
    PASS();
}

/* ---- L5: Fast Transforms ---- */
static void test_fast_transforms(void) {
    TEST("Walsh-Hadamard transform constant function");
    double f[4] = {1.0, 1.0, 1.0, 1.0};
    walsh_hadamard_transform(f, 2);
    /* f_hat(0) = 4, others 0 */
    assert(fabs(f[0] - 4.0) < 1e-9);
    assert(fabs(f[1]) < 1e-9);
    assert(fabs(f[2]) < 1e-9);
    assert(fabs(f[3]) < 1e-9);
    PASS();

    TEST("Walsh-Hadamard transform delta function");
    double g[4] = {1.0, 0.0, 0.0, 0.0}; /* delta at 0 */
    walsh_hadamard_transform(g, 2);
    /* All coefficients should be 1 */
    for (int32_t i = 0; i < 4; i++)
        assert(fabs(g[i] - 1.0) < 1e-9);
    PASS();

    TEST("poly_evaluate_all_bool");
    polynomial_t *p = poly_create(2, 10, false);
    poly_add_term(p, 1.0, 0x0);  /* constant 1 */
    poly_add_term(p, 2.0, 0x1);  /* 2*x0 */
    /* p = 1 + 2*x0, evaluations: [1, 3, 1, 3] */
    double evals[4];
    poly_evaluate_all_bool(p, evals, 2);
    assert(fabs(evals[0] - 1.0) < 1e-9);  /* x=(0,0) */
    assert(fabs(evals[1] - 3.0) < 1e-9);  /* x=(1,0) */
    assert(fabs(evals[2] - 1.0) < 1e-9);  /* x=(0,1) */
    assert(fabs(evals[3] - 3.0) < 1e-9);  /* x=(1,1) */
    poly_destroy(p);
    PASS();

    TEST("poly_and subset convolution");
    /* f = x0, g = x1 -> f AND g = 0 (never both true at same var?) */
    /* Actually f AND g = x0*x1 (pointwise AND) */
    polynomial_t *fa = poly_create(2, 10, true);
    polynomial_t *ga = poly_create(2, 10, true);
    poly_add_term(fa, 1.0, 0x1);  /* x0 */
    poly_add_term(ga, 1.0, 0x2);  /* x1 */
    polynomial_t *ha = poly_and(fa, ga);
    assert(ha != NULL);
    /* x0 AND x1 = x0*x1: the monomial for var mask 0x3 should exist */
    bool found_x0x1 = false;
    for (int32_t i = 0; i < ha->num_terms; i++) {
        if (ha->terms[i].mono.vars == 0x3 && fabs(ha->terms[i].coeff) > 0.5)
            found_x0x1 = true;
    }
    assert(found_x0x1);
    poly_destroy(fa); poly_destroy(ga); poly_destroy(ha);
    PASS();
}

/* ---- L6: OV Problem Solvers ---- */
static void test_ov_solvers(void) {
    TEST("ov_solve_brute_force basic");
    ov_instance_t *A = ov_create(3, 4);
    ov_instance_t *B = ov_create(3, 4);
    ov_set_vector(A, 0, 0x3);  /* 0011 */
    ov_set_vector(A, 1, 0xC);  /* 1100 */
    ov_set_vector(A, 2, 0xF);  /* 1111 */
    ov_set_vector(B, 0, 0xC);  /* 1100 -> orthogonal to A[0] */
    ov_set_vector(B, 1, 0x1);  /* 0001 -> dot A[1]=1100 = 0? No: 0001 & 1100 = 0 => orthogonal! */
    ov_set_vector(B, 2, 0x0);  /* 0000 -> orthogonal to everything */
    ov_result_t r = ov_solve_brute_force(A, B);
    assert(r.orthogonal_pair_exists);
    ov_destroy(A); ov_destroy(B);
    PASS();

    TEST("ov_solve_packed same as brute force");
    ov_instance_t *A2 = ov_create(5, 8);
    ov_instance_t *B2 = ov_create(5, 8);
    for (int32_t i = 0; i < 5; i++) {
        ov_set_vector(A2, i, (uint64_t)(1 << (i % 8)));
        ov_set_vector(B2, i, (uint64_t)(~(1 << (i % 8)) & 0xFF));
    }
    ov_result_t r1 = ov_solve_brute_force(A2, B2);
    ov_result_t r2 = ov_solve_packed(A2, B2);
    assert(r1.orthogonal_pair_exists == r2.orthogonal_pair_exists);
    ov_destroy(A2); ov_destroy(B2);
    PASS();

    TEST("ov_count_orthogonal_pairs");
    ov_instance_t *A3 = ov_create(2, 2);
    ov_instance_t *B3 = ov_create(2, 2);
    ov_set_vector(A3, 0, 0x1); /* 01 */
    ov_set_vector(A3, 1, 0x2); /* 10 */
    ov_set_vector(B3, 0, 0x2); /* 10 -> orthogonal to A[0], not A[1] */
    ov_set_vector(B3, 1, 0x1); /* 01 -> orthogonal to A[1], not A[0] */
    int64_t cnt = ov_count_orthogonal_pairs(A3, B3);
    assert(cnt == 2);
    ov_destroy(A3); ov_destroy(B3);
    PASS();

    TEST("ov_verify_orthogonal");
    ov_instance_t *A4 = ov_create(3, 4);
    ov_instance_t *B4 = ov_create(3, 4);
    ov_set_vector(A4, 0, 0x3); ov_set_vector(B4, 0, 0xC);
    assert(ov_verify_orthogonal(A4, B4, 0, 0));  /* 0011 and 1100 */
    ov_set_vector(A4, 1, 0xF); ov_set_vector(B4, 1, 0xF);
    assert(!ov_verify_orthogonal(A4, B4, 1, 1)); /* 1111 and 1111 */
    ov_destroy(A4); ov_destroy(B4);
    PASS();
}

/* ---- L6: Williams SAT Algorithm ---- */
static void test_williams_sat(void) {
    TEST("williams_cnf_evaluate");
    williams_cnf_t *f = williams_cnf_create(3, 1, 3);
    int32_t clause[] = {1, 2, 3};  /* x1 OR x2 OR x3 */
    williams_cnf_add_clause(f, clause, 3);
    williams_assign_t *a = williams_assign_create(3);
    a->values[0] = 0; a->values[1] = 0; a->values[2] = 1;
    assert(williams_cnf_evaluate(f, a));
    a->values[0] = 0; a->values[1] = 0; a->values[2] = 0;
    assert(!williams_cnf_evaluate(f, a));
    williams_assign_destroy(a);
    williams_cnf_destroy(f);
    PASS();

    TEST("williams_brute_force satisfiable");
    williams_cnf_t *f2 = williams_cnf_create(4, 3, 3);
    int32_t c1[] = {1, 2, 3};
    int32_t c2[] = {-1, -2, 4};
    int32_t c3[] = {1, -3, -4};
    williams_cnf_add_clause(f2, c1, 3);
    williams_cnf_add_clause(f2, c2, 3);
    williams_cnf_add_clause(f2, c3, 3);
    williams_result_t r = williams_brute_force(f2);
    assert(r.satisfiable);
    assert(williams_verify_witness(f2, r.witness));
    williams_cnf_destroy(f2);
    free(r.witness);
    PASS();

    TEST("williams_solve matches brute force");
    williams_cnf_t *f3 = williams_cnf_create(6, 3, 3);
    int32_t d1[] = {1, 2, 3};
    int32_t d2[] = {-2, -4, 5};
    int32_t d3[] = {3, -5, 6};
    williams_cnf_add_clause(f3, d1, 3);
    williams_cnf_add_clause(f3, d2, 3);
    williams_cnf_add_clause(f3, d3, 3);
    williams_result_t rb = williams_brute_force(f3);
    williams_result_t rw = williams_solve(f3);
    assert(rb.satisfiable == rw.satisfiable);
    if (rb.satisfiable) {
        assert(williams_verify_witness(f3, rb.witness));
        assert(williams_verify_witness(f3, rw.witness));
    }
    williams_cnf_destroy(f3);
    free(rb.witness); free(rw.witness);
    PASS();
}

/* ---- L4/L8: Lower Bounds via Polynomial Degree ---- */
static void test_lower_bounds(void) {
    TEST("MOD_3 degree over GF(2)");
    uint8_t *tt_mod3 = build_mod_m_truth_table(4, 3);
    int32_t deg = gf_degree_from_truth_table(tt_mod3, 4, 2);
    /* MOD_3 over GF(2) should have high degree */
    assert(deg >= 1);
    free(tt_mod3);
    PASS();

    TEST("MOD_2 degree over GF(2) should be 1");
    uint8_t *tt_mod2 = build_mod_m_truth_table(4, 2);
    int32_t deg2 = gf_degree_from_truth_table(tt_mod2, 4, 2);
    assert(deg2 == 1);
    free(tt_mod2);
    PASS();

    TEST("Razborov degree bound computation");
    razborov_approx_params_t params = {3, 100, 4, 2};
    int32_t bound = razborov_degree_bound(&params, 0.1);
    assert(bound >= 1);
    PASS();

    TEST("Smolensky analysis MOD_3 vs GF(2)");
    smolensky_result_t sr = smolensky_analyze(6, 3, 2);
    assert(sr.requires_high_degree);  /* MOD_3 over GF(2) requires high degree */
    PASS();

    TEST("algebraic_circuit create/evaluate");
    algebraic_circuit_t *ac = ac_create(10, 3);
    int32_t x0 = ac_add_input(ac, 0);
    int32_t x1 = ac_add_input(ac, 1);
    int32_t sum = ac_add_add(ac, x0, x1);
    ac->output_gate = sum;
    double inputs[] = {1.0, 1.0, 0.0};
    double val = ac_evaluate(ac, inputs);
    assert(fabs(val - 2.0) < 1e-9);
    ac_destroy(ac);
    PASS();

    TEST("ac_compute_degree");
    algebraic_circuit_t *ac2 = ac_create(10, 2);
    int32_t xi0 = ac_add_input(ac2, 0);
    int32_t xi1 = ac_add_input(ac2, 1);
    int32_t prod = ac_add_mul(ac2, xi0, xi1);
    ac2->output_gate = prod;
    assert(ac_compute_degree(ac2) == 2);  /* x0*x1 has degree 2 */
    ac_destroy(ac2);
    PASS();

    TEST("circuit complexity estimate");
    uint8_t tt_and[8] = {0,0,0,0,0,0,0,1};
    circuit_complexity_estimate_t est = estimate_circuit_complexity(tt_and, 3);
    assert(est.exact_degree == 3);
    assert(est.sparsity >= 1);
    PASS();
}

int main(void) {
    printf("=== Polynomial Method Consequences Test Suite ===\n\n");

    printf("--- L1: Core Definitions ---\n");
    test_l1_definitions();

    printf("\n--- L3: GF(p) Arithmetic ---\n");
    test_gf_arithmetic();

    printf("\n--- L3: Polynomial Operations ---\n");
    test_polynomial_operations();

    printf("\n--- L4: Polynomial Representations ---\n");
    test_polynomial_representations();

    printf("\n--- L5: Fast Transforms ---\n");
    test_fast_transforms();

    printf("\n--- L6: OV Problem Solvers ---\n");
    test_ov_solvers();

    printf("\n--- L6: Williams SAT Algorithm ---\n");
    test_williams_sat();

    printf("\n--- L4/L8: Lower Bounds ---\n");
    test_lower_bounds();

    printf("\n==========================================\n");
    printf("  RESULTS: %d / %d tests passed\n", tests_passed, tests_total);
    printf("==========================================\n");

    return (tests_passed == tests_total) ? 0 : 1;
}
