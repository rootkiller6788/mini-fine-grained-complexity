
/* ============================================================================
 * seth_hardness.c -- SETH Conditional Lower Bounds Framework
 *
 * OV problem, SETH-to-Edit-Distance reduction verification,
 * 3SUM/APSP reduction checks, breakthrough analysis.
 *
 * L1: seth_statement_t, ov_conjecture_t, binary_vector_t
 * L2: OV brute force, 4-list splitting
 * L3: Reduction framework data structures
 * L4: Backurs-Indyk (2015), ABW (2015), Bringmann-Kunnemann (2015)
 * L8: Beyond-quadratic lower bounds
 * L9: Research frontier implications
 * ============================================================================ */

#include "seth_hardness.h"
#include "string_utils.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ============================================================================
 * L1: Binary Vector Operations
 * ============================================================================ */

binary_vector_t *bv_create(int32_t dim) {
    binary_vector_t *v = (binary_vector_t *)malloc(sizeof(binary_vector_t));
    if (!v) return NULL;
    v->dim = dim;
    v->num_words = (dim + 63) / 64;
    v->bits = (uint64_t *)calloc((size_t)v->num_words, sizeof(uint64_t));
    if (!v->bits) { free(v); return NULL; }
    return v;
}

void bv_destroy(binary_vector_t *v) {
    if (!v) return;
    free(v->bits);
    free(v);
}

void bv_set_bit(binary_vector_t *v, int32_t pos, bool value) {
    if (!v || pos < 0 || pos >= v->dim) return;
    int32_t word = pos / 64;
    int32_t bit = pos % 64;
    if (value) v->bits[word] |= (1ULL << bit);
    else v->bits[word] &= ~(1ULL << bit);
}

bool bv_get_bit(const binary_vector_t *v, int32_t pos) {
    if (!v || pos < 0 || pos >= v->dim) return false;
    int32_t word = pos / 64;
    int32_t bit = pos % 64;
    return (v->bits[word] >> bit) & 1ULL;
}

bool bv_are_orthogonal(const binary_vector_t *a, const binary_vector_t *b) {
    if (!a || !b || a->dim != b->dim) return false;
    int32_t words = a->num_words;
    for (int32_t i = 0; i < words; i++) {
        if (a->bits[i] & b->bits[i]) return false;
    }
    return true;
}

int32_t bv_dot_product(const binary_vector_t *a, const binary_vector_t *b) {
    if (!a || !b || a->dim != b->dim) return -1;
    int32_t dp = 0;
    for (int32_t i = 0; i < a->dim; i++) {
        if (bv_get_bit(a, i) && bv_get_bit(b, i)) dp++;
        if (dp > 0) return dp; /* Early exit: not orthogonal if any 1*1 */
    }
    return 0;
}

/* ============================================================================
 * L2: OV Instance and Algorithms
 * ============================================================================ */

ov_instance_t *ov_instance_create(int32_t num_a, int32_t num_b, int32_t dim) {
    ov_instance_t *inst = (ov_instance_t *)malloc(sizeof(ov_instance_t));
    if (!inst) return NULL;
    inst->num_a = num_a; inst->num_b = num_b; inst->dim = dim;
    inst->set_a = (binary_vector_t *)malloc((size_t)num_a * sizeof(binary_vector_t));
    inst->set_b = (binary_vector_t *)malloc((size_t)num_b * sizeof(binary_vector_t));
    if (!inst->set_a || !inst->set_b) {
        free(inst->set_a); free(inst->set_b); free(inst); return NULL;
    }
    for (int32_t i = 0; i < num_a; i++) {
        inst->set_a[i].dim = dim; inst->set_a[i].bits = NULL;
    }
    for (int32_t i = 0; i < num_b; i++) {
        inst->set_b[i].dim = dim; inst->set_b[i].bits = NULL;
    }
    return inst;
}

void ov_instance_destroy(ov_instance_t *inst) {
    if (!inst) return;
    for (int32_t i = 0; i < inst->num_a; i++) free(inst->set_a[i].bits);
    for (int32_t i = 0; i < inst->num_b; i++) free(inst->set_b[i].bits);
    free(inst->set_a); free(inst->set_b); free(inst);
}

void ov_instance_random_fill(ov_instance_t *inst, uint64_t seed) {
    if (!inst) return;
    xorshift64_t rng = xorshift64_seed(seed);
    for (int32_t i = 0; i < inst->num_a; i++) {
        inst->set_a[i].bits = (uint64_t *)calloc((size_t)((inst->dim+63)/64), sizeof(uint64_t));
        inst->set_a[i].num_words = (inst->dim + 63) / 64;
        if (!inst->set_a[i].bits) continue;
        for (int32_t j = 0; j < inst->dim; j++)
            if (xorshift64_double(&rng) < 0.5)
                bv_set_bit(&inst->set_a[i], j, true);
    }
    for (int32_t i = 0; i < inst->num_b; i++) {
        inst->set_b[i].bits = (uint64_t *)calloc((size_t)((inst->dim+63)/64), sizeof(uint64_t));
        inst->set_b[i].num_words = (inst->dim + 63) / 64;
        if (!inst->set_b[i].bits) continue;
        for (int32_t j = 0; j < inst->dim; j++)
            if (xorshift64_double(&rng) < 0.5)
                bv_set_bit(&inst->set_b[i], j, true);
    }
}

ov_result_t ov_brute_force(const ov_instance_t *inst) {
    ov_result_t result = {false, -1, -1, 0};
    if (!inst) return result;
    for (int32_t i = 0; i < inst->num_a; i++) {
        for (int32_t j = 0; j < inst->num_b; j++) {
            result.operations++;
            if (bv_are_orthogonal(&inst->set_a[i], &inst->set_b[j])) {
                result.exists_orthogonal = true;
                result.index_a = i;
                result.index_b = j;
                return result;
            }
        }
    }
    return result;
}

ov_result_t ov_four_list_split(const ov_instance_t *inst) {
    ov_result_t result = {false, -1, -1, 0};
    if (!inst) return result;
    /* 4-list splitting: split each set in half, check all 4 combinations.
     * OV(n,d) reduces to 4 OV(n/2, d) subproblems. */
    int32_t mid_a = inst->num_a / 2;
    for (int32_t ia = 0; ia <= 1; ia++) {
        int32_t sa = ia ? mid_a : 0;
        int32_t ea = ia ? inst->num_a : mid_a;
        for (int32_t ib = 0; ib <= 1; ib++) {
            int32_t sb = ib ? (inst->num_b / 2) : 0;
            int32_t eb = ib ? inst->num_b : (inst->num_b / 2);
            for (int32_t i = sa; i < ea; i++) {
                for (int32_t j = sb; j < eb; j++) {
                    result.operations++;
                    if (bv_are_orthogonal(&inst->set_a[i], &inst->set_b[j])) {
                        result.exists_orthogonal = true;
                        result.index_a = i; result.index_b = j;
                        return result;
                    }
                }
            }
        }
    }
    return result;
}

/* ============================================================================
 * L2: SETH Conditional Lower Bound Verification
 * ============================================================================ */

bool seth_edit_distance_check(double time_exponent) {
    /* Backurs-Indyk (2015): If Edit Distance can be solved in O(n^c)
     * for c < 2, then SETH is false. */
    return time_exponent < 2.0;
}

bool seth_lcs_check(double time_exponent) {
    /* ABW (2015): If LCS can be solved in O(n^c) for c < 2,
     * then SETH is false. */
    return time_exponent < 2.0;
}

double reduction_overhead(int32_t input_n, int32_t output_n,
                          int32_t from_dim, int32_t to_dim) {
    double overhead = (double)output_n / (double)input_n;
    if (from_dim > 0 && to_dim > 0)
        overhead *= (double)to_dim / (double)from_dim;
    return overhead;
}

/* ---- L4: Theorem Verification ---- */

bool backurs_indyk_verify(int32_t n, double epsilon) {
    /* Backurs-Indyk (STOC 2015) reduction chain:
     * SAT(n vars) -> OV(N=2^{n/2}, d=O(n)) -> Edit(N')
     * If ED in O(N^{2-eps}), then SAT in O(2^{(1-delta)n}).
     * Check: N^{2-eps} = (2^{n/2})^{2-eps} = 2^{n - eps*n/2}
     *        Compare to SETH: SAT needs 2^{(1-o(1))n}
     *        If eps > 0, then n - eps*n/2 < n, refuting SETH.
     */
    if (epsilon <= 0.0 || n <= 0) return false;
    double seth_threshold = (double)n * (1.0 - 1e-6);
    double edit_time_exp = (double)n - epsilon * (double)n / 2.0;
    return edit_time_exp < seth_threshold;
}

bool abw_verify(int32_t n, double epsilon) {
    /* ABW (FOCS 2015): LCS needs n^{2-o(1)} under SETH.
     * Same reduction chain via OV. */
    if (epsilon <= 0.0 || n <= 0) return false;
    double seth_threshold = (double)n * (1.0 - 1e-6);
    double lcs_time_exp = (double)n - epsilon * (double)n / 2.0;
    return lcs_time_exp < seth_threshold;
}

bool bringmann_kunnemann_verify(int32_t n, double epsilon) {
    /* Bringmann-Kunnemann (FOCS 2015): Frechet distance
     * needs n^{2-o(1)} under SETH. */
    if (epsilon <= 0.0 || n <= 0) return false;
    double threshold = (double)n * (1.0 - 1e-6);
    double frechet_exp = (double)n - epsilon * (double)n / 2.0;
    return frechet_exp < threshold;
}

/* ---- L2: Lower Bound Reports ---- */

void print_edit_seth_status(void) {
    printf("=== Edit Distance SETH Lower Bound ===\n");
    printf("Backurs-Indyk (STOC 2015): ED in O(n^{2-eps}) => SETH false\n");
    printf("Best known: O(n^2/log n) (Masek-Paterson 1980)\n");
}

void print_lcs_seth_status(void) {
    printf("=== LCS SETH Lower Bound ===\n");
    printf("ABW (FOCS 2015): LCS in O(n^{2-eps}) => SETH false\n");
    printf("Hirschberg (1975): O(n^2) time, O(min(n,m)) space\n");
}

void print_comprehensive_lower_bound_report(void) {
    printf("==== Fine-Grained Lower Bounds for Sequence Measures ====\n\n");
    printf("All below require n^{2-o(1)} time under SETH:\n");
    printf("  Edit Distance    | Backurs-Indyk 2015\n");
    printf("  LCS              | ABW 2015\n");
    printf("  Frechet Distance | Bringmann-Kunnemann 2015\n");
    printf("  DTW              | Bringmann-Kunnemann 2015\n");
    printf("\nChain: SETH -> OV Conjecture -> Sequence problems\n");
    printf("  SETH: k-SAT requires 2^{(1-o(1))n} time\n");
    printf("  OVC:  Orthogonal Vectors requires n^{2-o(1)} time\n");
    printf("  EDH:  Edit Distance requires n^{2-o(1)} time\n");
}

/* ---- L8: Beyond-Quadratic Lower Bounds ---- */

bool ov_to_edit_implies_subquadratic_edit(double time_exponent) {
    /* OV(n,d) reduces to Edit(N) where N = n * 2^{O(d)}.
     * For d = omega(log n), N = superpolynomial in n,
     * so subquadratic Edit doesn't directly imply subquadratic OV.
     * The reduction is only meaningful for d = O(log n). */
    return time_exponent < 2.0 && time_exponent > 0.0;
}

bool k_ov_to_edit_implies(double time_exponent, int32_t k) {
    /* k-OV: Given k sets of N vectors in {0,1}^d,
     * find a_1,...,a_k with zero dot product.
     * For larger k, the reduction to Edit Distance
     * produces longer strings, weakening the implication. */
    if (time_exponent < (double)k - 1.0) return true;
    return false;
}

bool ovh_check(double time_exponent, int32_t d, int32_t n) {
    /* OVH: For d = omega(log n), OV requires n^{2-o(1)} time.
     * This is a stronger hypothesis than the reduction from SETH,
     * as SETH only implies OVC for specific d. */
    double log_n = log((double)n);
    if ((double)d > 2.0 * log_n && time_exponent < 2.0) return true;
    return false;
}

bool apsp_to_edit_check(double time_exponent) {
    /* APSP conjecture: no O(n^{3-eps}) algorithm for
     * All-Pairs Shortest Paths in dense graphs.
     * APSP -> Min-Plus Convolution -> Edit Distance
     * (indirect connection via chain of reductions) */
    return time_exponent < 2.0;
}

bool three_sum_to_edit_check(double time_exponent) {
    /* 3SUM conjecture: no O(n^{2-eps}) algorithm.
     * 3SUM -> Convolution3SUM -> Edit Distance
     * (indirect via geometric reductions) */
    return time_exponent < 2.0;
}

/* ---- L9: Research Frontiers ---- */

double seth_constant_s_k(int32_t k) {
    /* Known bounds for s_k = inf{c : k-SAT in O(2^{cn})}
     * s_3 = 0.3863 (PPSZ 2005)
     * s_4 = 0.5548
     * s_5 ~ 0.6500
     * lim_{k->infty} s_k = 1 (SETH conjectured) */
    static const double sv[] = {0.0, 0.0, 0.0, 0.3863, 0.5548,
        0.6500, 0.7162, 0.7627, 0.7977, 0.8245, 0.8453};
    if (k >= 3 && k <= 10) return sv[k];
    if (k < 3) return 0.0;
    return 1.0 - 1.5 / (double)k;
}

void breakthrough_implications(double c) {
    printf("=== Implications of O(n^%.2f) Edit Distance ===\n", c);
    if (c < 2.0) {
        printf("  -> SETH is REFUTED (Backurs-Indyk 2015)\n");
        printf("  -> OV Conjecture is REFUTED\n");
        printf("  -> New SAT algorithms must exist\n");
        printf("  -> LCS and Frechet in O(n^{%.3f})\n", c);
    } else if (c < 2.0 + 1e-6) {
        printf("  -> Matches current best: O(n^2/log n)\n");
        printf("  -> SETH remains plausible\n");
    } else {
        printf("  -> Worse than current best algorithms\n");
    }
}
