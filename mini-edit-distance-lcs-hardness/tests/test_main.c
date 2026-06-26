/* ============================================================================
 * test_main.c -- Tests for Edit Distance, LCS, Alignment, Frechet, SETH
 * ============================================================================ */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "../include/edit_distance.h"
#include "../include/lcs.h"
#include "../include/seth_hardness.h"
#include "../include/alignment.h"
#include "../include/frechet_distance.h"
#include "../include/string_utils.h"

#define TEST(n) printf("  %-45s ", n)
#define PASS() do { printf("PASS\n"); tests_passed++; } while(0)
static int tests_passed = 0;

/* ---- L1-L2: Edit Distance ---- */
void test_ed_basic(void) {
    TEST("ED identical"); assert(edit_distance_full("abc","abc")==0); PASS();
    TEST("ED empty"); assert(edit_distance_full("","")==0); PASS();
    TEST("ED kitten->sitting=3"); assert(edit_distance_full("kitten","sitting")==3); PASS();
    TEST("ED saturday->sunday=3"); assert(edit_distance_full("Saturday","Sunday")==3); PASS();
    TEST("ED abc->xyz=3"); assert(edit_distance_full("abc","xyz")==3); PASS();
    TEST("ED single sub"); assert(edit_distance_full("a","b")==1); PASS();
    TEST("ED insert"); assert(edit_distance_full("a","ab")==1); PASS();
    TEST("ED delete"); assert(edit_distance_full("ab","a")==1); PASS();
}
void test_ed_weighted(void) {
    TEST("Weighted ED delete cheap");
    edit_cost_model_t c = {10,1,1,0};
    assert(edit_distance_weighted("abc","ab",&c)==1); PASS();
    TEST("Weighted ED sub expensive");
    c = (edit_cost_model_t){1,1,100,0};
    assert(edit_distance_weighted("a","b",&c)==2); PASS();
}
void test_ed_linear(void) {
    TEST("Linear space matches full");
    assert(edit_distance_linear_space("saturday","sunday")==edit_distance_full("saturday","sunday")); PASS();
}
void test_ed_bounded(void) {
    TEST("Bounded k=0 identical"); assert(edit_distance_bounded("abc","abc",0)==0); PASS();
    TEST("Bounded k=0 diff"); assert(edit_distance_bounded("abc","abd",0)==1); PASS();
}
void test_hamming(void) {
    TEST("Hamming abc vs abd=1"); assert(hamming_distance("abc","abd")==1); PASS();
    TEST("Hamming unequal len"); assert(hamming_distance("abc","ab")==-1); PASS();
}
void test_dl(void) {
    TEST("DL transposition abc->acb=1"); assert(damerau_levenshtein("abc","acb")==1); PASS();
}
void test_jw(void) {
    TEST("Jaro-Winkler similar"); assert(jaro_winkler("martha","marhta",0.1)>0.9); PASS();
    TEST("JW identical"); double j=jaro_winkler("abc","abc",0.1); assert(fabs(j-1.0)<1e-9); PASS();
}

/* ---- L3-L4: LCS ---- */
void test_lcs_len(void) {
    TEST("LCS AGGTAB,GXTXAYB=4"); assert(lcs_length("AGGTAB","GXTXAYB")==4); PASS();
    TEST("LCS identical"); assert(lcs_length("abc","abc")==3); PASS();
    TEST("LCS none"); assert(lcs_length("abc","xyz")==0); PASS();
}
void test_lcs_res(void) {
    TEST("LCS result GTAB");
    lcs_result_t *r=lcs_compute("AGGTAB","GXTXAYB");
    assert(r&&r->length==4&&strcmp(r->subsequence,"GTAB")==0);
    lcs_result_destroy(r); PASS();
}
void test_hirschberg(void) {
    TEST("Hirschberg matches standard");
    lcs_result_t *s=lcs_compute("AGGTAB","GXTXAYB"),*h=lcs_hirschberg("AGGTAB","GXTXAYB");
    assert(s->length==h->length&&strcmp(s->subsequence,h->subsequence)==0);
    lcs_result_destroy(s); lcs_result_destroy(h); PASS();
}
void test_lcs_var(void) {
    TEST("LCSubstring abcdef,zbcde=4"); assert(longest_common_substring("abcdef","zbcde")==4); PASS();
    TEST("SCS ab,ac=3"); assert(shortest_common_supersequence_length("ab","ac")==3); PASS();
}
void test_lis(void) {
    TEST("LIS works"); int32_t a[]={3,1,4,1,5,9,2,6,5,3,5};
    assert(longest_increasing_subsequence(a,11)==4); PASS();
}
void test_lcis(void) {
    TEST("LCIS works"); int32_t a[]={3,4,9,1},b[]={5,3,8,9,10,2,1};
    assert(longest_common_increasing_subsequence(a,4,b,7)==2); PASS();
}
void test_ed_via_lcs(void) {
    TEST("ED via LCS formula ok");
    int32_t evl = edit_via_lcs("abc", "abd");
    /* ED via LCS: |a|+|b|-2*LCS = 3+3-2*2 = 2 (subst=del+ins, cost=2 per sub) */
    assert(evl == 2); PASS();
}


/* ---- L5: Alignment ---- */
void test_nw(void) {
    TEST("NW alignment");
    alignment_scoring_t s={1,-1,-2,-1,false};
    alignment_result_t *r=align_needleman_wunsch("AGCT","AGGT",&s);
    assert(r&&r->num_columns>0); alignment_result_destroy(r); PASS();
}
void test_sw(void) {
    TEST("SW alignment"); alignment_scoring_t s={1,-1,-2,-1,false};
    alignment_result_t *r=align_smith_waterman("ACGT","ACGT",&s);
    assert(r); alignment_result_destroy(r); PASS();
}
void test_submat(void) {
    TEST("Substitution matrix identity");
    substitution_matrix_t *m=substitution_matrix_load("identity");
    assert(m&&substitution_score(m,'A','A')==1);
    substitution_matrix_destroy(m); PASS();
}
void test_bioinfo(void) {
    TEST("DNA rev comp"); char *rc=dna_reverse_complement("ATGC");
    assert(rc&&strcmp(rc,"GCAT")==0); free(rc); PASS();
    TEST("GC content"); double gc=dna_gc_content("ATGC");
    assert(fabs(gc-0.5)<1e-9); PASS();
    TEST("DNA to protein"); char *p=dna_to_protein("ATGGCCGGG");
    assert(p); free(p); PASS();
}

/* ---- L6: Frechet ---- */
void test_frechet(void) {
    TEST("Frechet identical"); double x[]={0,1,2},y[]={0,0,0};
    curve_t *c1=curve_from_arrays(x,y,3),*c2=curve_from_arrays(x,y,3);
    assert(fabs(discrete_frechet_distance(c1,c2))<1e-9);
    curve_destroy(c1); curve_destroy(c2); PASS();
    TEST("Frechet decision"); double qy2[]={0,0.5,0};
    c1=curve_from_arrays(x,y,3); c2=curve_from_arrays(x,qy2,3);
    assert(frechet_decision(c1,c2,5.0));
    curve_destroy(c1); curve_destroy(c2); PASS();
    TEST("DTW finite"); double py[]={0,1,0};
    c1=curve_from_arrays(x,py,3); c2=curve_from_arrays(x,qy2,3);
    assert(isfinite(dtw_distance(c1,c2)));
    curve_destroy(c1); curve_destroy(c2); PASS();
}

/* ---- L7: SETH ---- */
void test_seth(void) {
    TEST("SETH ED check c=1.5"); assert(seth_edit_distance_check(1.5)); PASS();
    TEST("SETH LCS check c=1.9"); assert(seth_lcs_check(1.9)); PASS();
    TEST("Backurs-Indyk"); assert(backurs_indyk_verify(100,0.1)); PASS();
    TEST("ABW"); assert(abw_verify(100,0.1)); PASS();
    TEST("Bringmann-Kunnemann"); assert(bringmann_kunnemann_verify(100,0.1)); PASS();
}
void test_ov(void) {
    TEST("OV vector create"); binary_vector_t *v=bv_create(128); assert(v&&v->dim==128); bv_destroy(v); PASS();
    TEST("OV bits"); v=bv_create(64); bv_set_bit(v,0,1); assert(bv_get_bit(v,0)); bv_destroy(v); PASS();
    TEST("OV orthogonal"); binary_vector_t *a=bv_create(64),*b2=bv_create(64);
    bv_set_bit(a,0,1); bv_set_bit(b2,1,1); assert(bv_are_orthogonal(a,b2));
    bv_set_bit(b2,0,1); assert(!bv_are_orthogonal(a,b2));
    bv_destroy(a); bv_destroy(b2); PASS();
    TEST("OV brute"); ov_instance_t *inst=ov_instance_create(10,10,8);
    ov_instance_random_fill(inst,12345); ov_result_t res=ov_brute_force(inst);
    assert(res.operations>0); ov_instance_destroy(inst); PASS();
}
void test_seth_const(void) {
    TEST("s_k monotonic"); double s3=seth_constant_s_k(3);
    assert(fabs(s3-0.3863)<0.01); assert(seth_constant_s_k(10)>s3); PASS();
}


/* ---- L8: String Utils & Advanced ---- */
void test_str_utils(void) {
    TEST("str_len"); assert(str_len("hello")==5); PASS();
    TEST("str_eq"); assert(str_eq("abc","abc")&&!str_eq("abc","ABC")); PASS();
    TEST("str_reverse"); char *s=str_dup("hello"); str_reverse(s);
    assert(strcmp(s,"olleh")==0); free(s); PASS();
    TEST("palindrome"); assert(str_is_palindrome("racecar")&&!str_is_palindrome("hello")); PASS();
    TEST("random DNA"); char *r=dna_random(20,42); assert(r&&str_len(r)==20); free(r); PASS();
    TEST("str_mutate"); char *m=str_mutate("hello",0.1,42); assert(m); free(m); PASS();
}
void test_fr(void) {
    TEST("Four Russians ED"); assert(edit_distance_four_russians("abcdefgh","abxyefgh")>=0); PASS();
    TEST("Four Russians LCS"); assert(lcs_four_russians("abcdefghijkl","abxyefghpqrs")>=0); PASS();
}

/* ---- Main test runner ---- */
int main(void) {
    printf("=== mini-edit-distance-lcs-hardness: Test Suite ===\n\n");
    printf("[L1-L2] Edit Distance\n"); test_ed_basic(); test_ed_weighted();
    test_ed_linear(); test_ed_bounded(); test_hamming(); test_dl(); test_jw();
    printf("\n[L3-L4] LCS\n"); test_lcs_len(); test_lcs_res();
    test_hirschberg(); test_lcs_var(); test_lis(); test_lcis(); test_ed_via_lcs();
    printf("\n[L5] Alignment\n"); test_nw(); test_sw(); test_submat(); test_bioinfo();
    printf("\n[L6] Frechet\n"); test_frechet();
    printf("\n[L7] SETH Hardness\n"); test_seth(); test_ov(); test_seth_const();
    printf("\n[L8] Utilities + Advanced\n"); test_str_utils(); test_fr();
    printf("\n========================================\n");
    printf("All %d tests PASSED\n", tests_passed);
    return 0;
}

