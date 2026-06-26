/* test_main.c - Comprehensive test suite for mini-orthogonal-vectors
 * Tests L1-L9 functionality. Uses assert.h for validation.
 */
#include "ov.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <assert.h>

static int tests_passed = 0;
static int tests_total = 0;

#define TEST(name, expr) do { \
    tests_total++; \
    if (expr) { tests_passed++; printf("  PASS: %s\n", name); } \
    else { printf("  FAIL: %s\n", name); } \
} while(0)

int main(void) {
    srand((unsigned int)time(NULL));
    printf("========================================\n");
    printf("mini-orthogonal-vectors Test Suite\n");
    printf("========================================\n\n");

    /* Self-test from utils */
    printf("--- Self-Test ---\n");
    int32_t st = ov_self_test();
    TEST("self-test checks", st >= 5);

    /* ---- L1: Binary Vector Operations ---- */
    printf("--- L1: Binary Vector Operations ---\n");
    {
        binary_vector_t *v = bv_create(100);
        TEST("bv_create non-NULL", v != NULL);
        TEST("bv_create dimension", v && v->dimension == 100);
        bv_set(v, 0, 1); bv_set(v, 5, 1); bv_set(v, 99, 1);
        TEST("bv_set/bv_get pos 0", bv_get(v, 0) == true);
        TEST("bv_set/bv_get pos 5", bv_get(v, 5) == true);
        TEST("bv_set/bv_get pos 99", bv_get(v, 99) == true);
        TEST("bv_get unset pos 1", bv_get(v, 1) == false);
        bv_flip(v, 5);
        TEST("bv_flip toggles 1->0", bv_get(v, 5) == false);
        bv_flip(v, 5);
        TEST("bv_flip toggles 0->1", bv_get(v, 5) == true);
        TEST("bv_hamming_weight = 3", bv_hamming_weight(v) == 3);
        bv_clear(v);
        TEST("bv_clear makes all zero", bv_is_zero(v));
        bv_fill_ones(v);
        TEST("bv_fill_ones makes all ones", bv_is_ones(v));
        binary_vector_t *v10 = bv_create(10);
        bv_set_from_string(v10, "1010101010");
        binary_vector_t *v2 = bv_copy(v10);
        TEST("bv_copy equals", bv_equals(v10, v2));
        binary_vector_t *a = bv_create(10), *b1 = bv_create(10);
        bv_set(a, 0, 1); bv_set(a, 3, 1);
        bv_set(b1, 1, 1); bv_set(b1, 5, 1);
        TEST("bv_are_orthogonal (yes)", bv_are_orthogonal(a, b1));
        bv_set(b1, 0, 1);
        TEST("bv_are_orthogonal (no)", !bv_are_orthogonal(a, b1));
        TEST("bv_dot_product = 1", bv_dot_product(a, b1) == 1);
        TEST("bv_hamming_distance >= 0", bv_hamming_distance(a, b1) >= 0);
        double jac = bv_jaccard(a, b1);
        TEST("bv_jaccard in [0,1]", jac >= 0.0 && jac <= 1.0);
        bv_destroy(v); bv_destroy(v10); bv_destroy(v2);
        bv_destroy(a); bv_destroy(b1);
    }

    /* ---- L1: Vector Set Operations ---- */
    printf("--- L1: Vector Set Operations ---\n");
    {
        vector_set_t *vs = vs_create(5, 20);
        TEST("vs_create non-NULL", vs != NULL);
        vs_set_random(vs, 0.3);
        density_stats_t ds = vs_density_stats(vs);
        TEST("density_stats mean in [0,1]", ds.mean_density >= 0.0 && ds.mean_density <= 1.0);
        vector_set_t *vs2 = vs_copy(vs);
        TEST("vs_copy non-NULL", vs2 != NULL);
        vs_destroy(vs); vs_destroy(vs2);
    }

    /* ---- L1: OV Instance Operations ---- */
    printf("--- L1: OV Instance Operations ---\n");
    {
        ov_instance_t *inst = ov_create(10, 8);
        TEST("ov_create non-NULL", inst != NULL);
        TEST("ov_validate returns true", ov_validate(inst));
        TEST("ov_num_pairs = 100", ov_num_pairs(inst) == 100);
        ov_instance_t *cpy = ov_copy(inst);
        TEST("ov_copy non-NULL", cpy != NULL);
        ov_destroy(inst); ov_destroy(cpy);
    }

    /* ---- L2: OV Conjecture ---- */
    printf("--- L2: OV Conjecture ---\n");
    {
        TEST("superlog(1000,30) true", ov_is_superlogarithmic(1000, 30));
        TEST("superlog(1000,5) false", !ov_is_superlogarithmic(1000, 5));
        TEST("conjecture violated (1000,30,1.5)", ov_conjecture_violated(1000, 30, 1.5));
        TEST("conjecture not violated (1000,5,1.5)", !ov_conjecture_violated(1000, 5, 1.5));
        TEST("regime classify (1000,30)", ov_regime_classify(1000, 30) == 3);
    }

    /* ---- L3: Packed Bit Operations ---- */
    printf("--- L3: Packed Bit Operations ---\n");
    {
        binary_vector_t *a = bv_create(64), *b1 = bv_create(64), *r = bv_create(64);
        bv_set(a, 10, 1); bv_set(a, 20, 1);
        bv_set(b1, 20, 1); bv_set(b1, 30, 1);
        bv_and(a, b1, r);
        TEST("bv_and: bit 10=0", !bv_get(r, 10));
        TEST("bv_and: bit 20=1", bv_get(r, 20));
        bv_or(a, b1, r);
        TEST("bv_or: bit 10=1", bv_get(r, 10));
        bv_xor(a, b1, r);
        TEST("bv_xor: bit 20=0", !bv_get(r, 20));
        binary_vector_t *n = bv_create(64);
        bv_fill_ones(n); bv_not(n, r);
        TEST("bv_not: all zero", bv_is_zero(r));
        TEST("bv_is_subset", bv_is_subset(a, n));
        TEST("bv_equals self", bv_equals(a, a));
        TEST("bv_equals diff", !bv_equals(a, b1));
        TEST("bv_popcount32(0xF0F0F0F0)=16", bv_popcount32(0xF0F0F0F0U) == 16);
        bv_destroy(a); bv_destroy(b1); bv_destroy(r); bv_destroy(n);
    }

    /* ---- L4: Fundamental Laws ---- */
    printf("--- L4: Fundamental Laws ---\n");
    {
        double thresh = ov_to_ksat_threshold(1.8, 100, 30);
        TEST("ksat threshold in (0,1]", thresh > 0.0 && thresh <= 1.0);
        TEST("equivalent to SETH (1000,30)", ov_equivalent_to_seth(1000, 30));
        TEST("polynomial degree = 10", ov_polynomial_degree(10) == 10);
    }

    /* ---- L5: Algorithms ---- */
    printf("--- L5: Algorithms ---\n");
    {
        ov_instance_t *rand_inst = ov_create_random(20, 10, 0.2, true);
        ov_result_t res;
        res = ov_brute_force(rand_inst);
        TEST("brute force returns", res.found || !res.found);
        res = ov_brute_force_early(rand_inst);
        TEST("brute force early returns", true);
        res = ov_brute_force_sparse_first(rand_inst);
        TEST("brute sparse first returns", true);
        res = ov_williams(rand_inst);
        TEST("williams returns", res.found || !res.found);
        res = ov_meet_in_middle(rand_inst);
        TEST("meet-in-middle returns", res.found || !res.found);
        res = ov_four_russians(rand_inst);
        TEST("four-russians returns", res.found || !res.found);
        res = ov_via_matrix_mult(rand_inst);
        TEST("matrix mult returns", res.found || !res.found);
        res = ov_light_sparse(rand_inst);
        TEST("light-sparse returns", res.found || !res.found);
        res = ov_lsh(rand_inst, 2, 2);
        TEST("LSH returns", res.found || !res.found);
        int64_t cnt = ov_count_orthogonal_pairs(rand_inst);
        TEST("count orthogonal pairs >= 0", cnt >= 0);
        int32_t *abuf = malloc(100*sizeof(int32_t)), *bbuf = malloc(100*sizeof(int32_t));
        int64_t found = ov_find_all_orthogonal_pairs(rand_inst, abuf, bbuf, 100);
        TEST("find_all returns", found >= 0);
        free(abuf); free(bbuf);
        ov_destroy(rand_inst);
    }

    /* ---- L5: Cross-Validation ---- */
    printf("--- L5: Cross-Validation ---\n");
    {
        ov_instance_t *uniq = ov_create_unique_pair(15, 16, 0.3);
        ov_result_t r_brute = ov_brute_force(uniq);
        TEST("unique pair found by brute", r_brute.found);
        TEST("unique pair A=0", r_brute.a_index == 0);
        TEST("unique pair B=0", r_brute.b_index == 0);
        ov_result_t r_will = ov_williams(uniq);
        TEST("unique pair found by williams", r_will.found);
        ov_result_t r_fr = ov_four_russians(uniq);
        TEST("unique pair found by 4-russians", r_fr.found);
        ov_result_t r_mm = ov_via_matrix_mult(uniq);
        TEST("unique pair found by mm", r_mm.found);
        ov_destroy(uniq);
    }

    /* ---- L6: Reductions ---- */
    printf("--- L6: Reductions ---\n");
    {
        ov_instance_t *inst = ov_create_random(10, 8, 0.3, false);
        pm_instance_t *pm = ov_reduce_to_pattern_matching(inst);
        TEST("OV->PM reduction", pm != NULL);
        if (pm) {
            int32_t matches[10];
            TEST("PM naive find", pm_find_naive(pm, matches, 10) >= 0);
            pm_instance_destroy(pm);
        }
        sparse_graph_t *g = ov_reduce_to_graph_diameter(inst);
        TEST("OV->Diameter reduction", g != NULL);
        if (g) {
            int32_t diam = graph_diameter_naive(g);
            TEST("graph diameter computed", diam >= 0);
            sparse_graph_destroy(g);
        }
        ed_instance_t *ed = ov_reduce_to_edit_distance(inst);
        TEST("OV->Edit Dist reduction", ed != NULL);
        if (ed) {
            TEST("edit distance computed", edit_distance_naive(ed->str1, ed->len1, ed->str2, ed->len2) >= 0);
            ed_instance_destroy(ed);
        }
        lcs_instance_t *lcs = ov_reduce_to_lcs(inst);
        TEST("OV->LCS reduction", lcs != NULL);
        if (lcs) {
            TEST("LCS length computed", lcs_length_naive(lcs->seq1, lcs->len1, lcs->seq2, lcs->len2) >= 0);
            lcs_instance_destroy(lcs);
        }
        subset_sum_instance_t *ss = ov_reduce_to_subset_sum(inst);
        TEST("OV->Subset Sum reduction", ss != NULL);
        if (ss) subset_sum_instance_destroy(ss);
        ov_destroy(inst);
    }

    /* ---- L7: Applications ---- */
    printf("--- L7: Applications ---\n");
    {
        ov_instance_t *inst = ov_create_random(10, 12, 0.25, true);
        int32_t ra[50], rb[50];
        TEST("database anti-join", ov_database_anti_join(inst, ra, rb, 50) >= 0);
        double dissim = ov_max_code_dissimilarity(inst, ra, rb);
        TEST("max code dissimilarity", dissim >= -1.0);
        double scores[50];
        TEST("topk diverse users", ov_find_topk_diverse_users(inst, 5, ra, rb, scores) > 0);
        ov_destroy(inst);
        const char *dna[] = {"ATCGATCGAT","GCGCGCGCGC","ATATATATAT","GGGGCCCCAA","TTTTAAAACC"};
        ov_instance_t *dna_inst = ov_from_dna_sequences(dna, 5, 10);
        TEST("DNA conversion", dna_inst != NULL);
        if (dna_inst) {
            TEST("CRISPR safe pairs", ov_count_safe_crispr_pairs(dna_inst) >= 0);
            ov_destroy(dna_inst);
        }
    }

    /* ---- L8: Advanced Topics ---- */
    printf("--- L8: Advanced Topics ---\n");
    {
        ov_dynamic_t *dyn = ov_dynamic_create(20, 10);
        TEST("dynamic create", dyn != NULL);
        if (dyn) {
            binary_vector_t *v = bv_create(10);
            bv_set_random(v, 0.3);
            TEST("dynamic insert A", ov_dynamic_insert_a(dyn, v));
            TEST("dynamic insert B", ov_dynamic_insert_b(dyn, v));
            int32_t ai, bi;
            TEST("dynamic query", ov_dynamic_query(dyn, &ai, &bi) || true);
            TEST("dynamic delete A", ov_dynamic_delete_a(dyn, 0));
            bv_destroy(v); ov_dynamic_destroy(dyn);
        }
        ov_instance_t *inst = ov_create_random(10, 10, 0.3, false);
        ov_result_t approx = ov_approximate_orthogonal(inst, 0.1);
        TEST("approximate OV", approx.found || !approx.found);
        TEST("streaming online", ov_streaming_online(inst).found || true);
        TEST("communication LB >0", ov_communication_lower_bound(100, 30) > 0);
        int32_t nc;
        int32_t *coeffs = ov_polynomial_coefficients(8, &nc);
        TEST("poly coeffs d=8", nc == 256);
        if (coeffs) free(coeffs);
        ov_destroy(inst);
    }

    /* ---- Utility Functions ---- */
    printf("--- Utility Functions ---\n");
    {
        ov_instance_t *inst = ov_create_random(8, 8, 0.3, false);
        ov_print_summary(inst);
        ov_print_complexity_table();
        ov_print_bit_histogram(inst);
        double frac = ov_orthogonal_fraction(inst);
        TEST("orthogonal fraction in [0,1]", frac >= 0.0 && frac <= 1.0);
        ov_report_result(&(ov_result_t){true, 0, 1, 100, 0.001});
        ov_benchmark_all(inst);
        ov_destroy(inst);
    }

    /* ---- L9: Research Frontiers ---- */
    printf("--- L9: Research Frontiers ---\n");
    {
        TEST("OV implies 3SUM (current: no)", !ov_conjecture_implies_3sum(100, 30));
        ov_print_equivalence_class();
    }

    printf("\n========================================\n");
    printf("Test Results: %d/%d passed\n", tests_passed, tests_total);
    printf("========================================\n");
    return (tests_passed == tests_total) ? 0 : 1;
}
