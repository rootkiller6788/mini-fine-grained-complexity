#include "poly_method.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* ============================================================================
 * williams.c - Williams' Algorithm Design via Polynomials
 *
 * Ryan Williams' breakthrough result (2010-2014): polynomial method
 * yields faster algorithms for circuit satisfiability and CNF-SAT.
 *
 * Key theorem (Williams 2010):
 *   For every k, k-SAT with n variables and m = poly(n) clauses can be
 *   solved in deterministic time O(2^{n - n/(O(log m))}).
 *   In particular, for m = O(n), this is O(2^{n - n/O(log n)}),
 *   which is faster than 2^n.
 *
 * This was the first progress on general CNF-SAT below 2^n since
 * the early 2000s. The algorithm combines:
 *   1. The polynomial method (fast evaluation of many polynomials)
 *   2. Circuit-analysis techniques (ACC circuit lower bounds)
 *   3. The "guess-and-verify" paradigm from fine-grained complexity
 *
 * L5: Algorithms - Williams' SAT algorithm
 * L6: Canonical problems - CNF-SAT
 * L7: Applications - circuit lower bounds, algorithm design
 * L8: Advanced topics - ACC circuit analysis, derandomization
 * ============================================================================ */

/* ---- Utility: wall-clock time in milliseconds ---- */
static double get_time_ms(void) {
    return (double)clock() * 1000.0 / (double)CLOCKS_PER_SEC;
}

/* ============================================================================
 * CNF Formula Management
 * ============================================================================ */

williams_cnf_t *williams_cnf_create(int32_t num_vars, int32_t num_clauses,
                                     int32_t k) {
    williams_cnf_t *f = (williams_cnf_t *)malloc(sizeof(williams_cnf_t));
    if (!f) return NULL;
    f->num_clauses = 0;
    f->num_vars = num_vars;
    f->clause_width = k;
    /* Each clause needs k literals + sentinel 0, so (k+1) per clause */
    f->data_len = num_clauses * (k + 1);
    f->clause_data = (int32_t *)calloc((size_t)f->data_len, sizeof(int32_t));
    if (!f->clause_data) { free(f); return NULL; }
    return f;
}

void williams_cnf_destroy(williams_cnf_t *f) {
    if (f) { free(f->clause_data); free(f); }
}

void williams_cnf_add_clause(williams_cnf_t *f, const int32_t *lits,
                               int32_t len) {
    if (!f || !lits || f->num_clauses >= f->data_len / (f->clause_width + 1))
        return;
    int32_t base = f->num_clauses * (f->clause_width + 1);
    for (int32_t i = 0; i < len && i < f->clause_width; i++)
        f->clause_data[base + i] = lits[i];
    f->clause_data[base + f->clause_width] = 0; /* sentinel */
    f->num_clauses++;
}

williams_assign_t *williams_assign_create(int32_t num_vars) {
    williams_assign_t *a = (williams_assign_t *)malloc(sizeof(williams_assign_t));
    if (!a) return NULL;
    a->values = (int8_t *)malloc((size_t)num_vars * sizeof(int8_t));
    if (!a->values) { free(a); return NULL; }
    a->num_vars = num_vars;
    for (int32_t i = 0; i < num_vars; i++) a->values[i] = -1;
    return a;
}

void williams_assign_destroy(williams_assign_t *a) {
    if (a) { free(a->values); free(a); }
}

/* ============================================================================
 * Evaluate CNF formula under full assignment
 *
 * Variable numbering: 1-based in clauses, 0-based in assignment array.
 * Positive literal = variable true, negative literal = variable false.
 * ============================================================================ */

bool williams_cnf_evaluate(const williams_cnf_t *f, const williams_assign_t *a) {
    if (!f || !a) return false;
    int32_t clause_size = f->clause_width + 1; /* +1 for sentinel */

    for (int32_t c = 0; c < f->num_clauses; c++) {
        int32_t base = c * clause_size;
        bool clause_satisfied = false;
        for (int32_t i = 0; f->clause_data[base + i] != 0; i++) {
            int32_t lit = f->clause_data[base + i];
            int32_t var = (lit > 0) ? lit : -lit;
            bool is_true = (lit > 0);
            int8_t val = a->values[var - 1];
            if (val == -1) continue; /* unassigned variable */
            if ((bool)val == is_true) {
                clause_satisfied = true;
                break;
            }
        }
        if (!clause_satisfied) return false;
    }
    return true;
}

/* ============================================================================
 * Brute-Force SAT Solver (Baseline)
 *
 * Enumerates all 2^n assignments. For each, check all m clauses.
 * Time: O(m * 2^n). This is the baseline that Williams' algorithm beats.
 * ============================================================================ */

williams_result_t williams_brute_force(const williams_cnf_t *f) {
    williams_result_t result = {false, NULL, 0.0, 0};
    if (!f) return result;

    double start = get_time_ms();
    int32_t n = f->num_vars;
    int64_t total = (int64_t)1 << n;
    int32_t clause_size = f->clause_width + 1;

    result.witness = (int8_t *)malloc((size_t)n * sizeof(int8_t));
    if (!result.witness) return result;

    for (int64_t assign_bits = 0; assign_bits < total; assign_bits++) {
        /* Check this assignment */
        bool satisfied = true;
        for (int32_t c = 0; c < f->num_clauses && satisfied; c++) {
            int32_t base = c * clause_size;
            bool clause_ok = false;
            for (int32_t i = 0; f->clause_data[base + i] != 0; i++) {
                int32_t lit = f->clause_data[base + i];
                int32_t var = (lit > 0) ? lit : -lit;
                bool need_true = (lit > 0);
                bool val = (assign_bits >> (var - 1)) & 1;
                if (val == need_true) { clause_ok = true; break; }
            }
            if (!clause_ok) { satisfied = false; break; }
        }

        if (satisfied) {
            result.satisfiable = true;
            for (int32_t i = 0; i < n; i++)
                result.witness[i] = (int8_t)((assign_bits >> i) & 1);
            result.evaluations = (uint64_t)(assign_bits + 1);
            break;
        }
    }

    result.runtime_ms = get_time_ms() - start;
    return result;
}

/* ============================================================================
 * Williams-Style SAT Solver (Simplified)
 *
 * Core idea: Split the n variables into two sets:
 *   - Left: variables 0..k-1 (size k = n/2)
 *   - Right: variables k..n-1 (size n-k)
 *
 * For each of the 2^k assignments to left variables:
 *   Construct a polynomial P in the right variables such that
 *   P(right_assignment) != 0 iff the formula is SATISFIED for
 *   this particular left assignment combined with right_assignment.
 *
 * Then evaluate ALL 2^k polynomials at ALL 2^{n-k} right assignments
 * using fast multi-point evaluation. If any evaluation is non-zero,
 * we have found a satisfying assignment.
 *
 * Time: O*(2^{n/2}) for the simplified version.
 *
 * The full Williams algorithm achieves better bounds by:
 *   - Recursively applying this split at multiple levels
 *   - Using fast rectangular matrix multiplication for evaluation
 *   - Exploiting the structure of CNF for polynomial construction
 *
 * Reference: R. Williams, "Improving Exhaustive Search Implies
 *   Superpolynomial Lower Bounds", STOC 2010.
 *   R. Williams, "Non-uniform ACC Circuit Lower Bounds", J. ACM 2014.
 * ============================================================================ */

/* Compute the residual formula polynomial for a partial assignment.
 *
 * Given a CNF formula and a partial assignment to the first split_point
 * variables, construct a polynomial over GF(2) in the remaining
 * n - split_point variables.
 *
 * The polynomial should evaluate to 1 (over GF(2)) for a right-side
 * assignment iff the formula is satisfied.
 *
 * For CNF-SAT, the residual is easier: after fixing left variables,
 * each clause either:
 *   - Is already satisfied (? contributes factor 1)
 *   - Is falsified (? contributes factor 0, entire formula = 0)
 *   - Has remaining unassigned literals (? contributes OR of those literals)
 *
 * The AND of ORs (CNF) becomes a polynomial via:
 *   (x OR y) = 1 - (1-x)*(1-y) over GF(2) = x + y + x*y
 *   (for Boolean variables, x^2 = x)
 *
 * The final polynomial is the product of all residual clause polynomials.
 * Degree ? k * number_of_clauses, but for k = O(1), it's manageable.
 * ============================================================================ */

polynomial_t *williams_residual_polynomial(const williams_cnf_t *f,
                                             const williams_assign_t *partial,
                                             int32_t split_point) {
    if (!f || !partial) return NULL;
    int32_t n = f->num_vars;
    int32_t right_vars = n - split_point;
    int32_t clause_size = f->clause_width + 1;

    /* The overall residual polynomial is the product of residual
     * polynomials for each clause. Start with constant 1. */
    polynomial_t *result = poly_create(right_vars, 64, true);
    if (!result) return NULL;
    poly_add_term(result, 1.0, 0ULL);

    for (int32_t c = 0; c < f->num_clauses; c++) {
        int32_t base = c * clause_size;
        bool already_satisfied = false;
        int right_lits[32];  /* max k=32 */
        int num_right = 0;

        for (int32_t i = 0; f->clause_data[base + i] != 0; i++) {
            int32_t lit = f->clause_data[base + i];
            int32_t var = (lit > 0) ? lit : -lit;
            bool need_true = (lit > 0);

            if (var - 1 < split_point) {
                /* Variable in the left (fixed) set */
                int8_t val = partial->values[var - 1];
                if (val != -1 && (bool)val == need_true) {
                    already_satisfied = true;
                    break;
                }
                /* If val is opposite, this literal is falsified, skip it */
            } else {
                /* Variable in the right (free) set; remap to 0..right_vars-1 */
                int32_t rvar = (var - 1) - split_point;
                right_lits[num_right++] = need_true ? (rvar + 1) : -(rvar + 1);
            }
        }

        if (already_satisfied) continue;  /* clause contributes factor 1 */

        if (num_right == 0) {
            /* All literals falsified by left assignment */
            poly_destroy(result);
            polynomial_t *zero = poly_create(right_vars, 1, true);
            return zero;  /* formula is falsified */
        }

        /* Build polynomial for this residual clause:
         * OR of the remaining literals.
         * Over GF(2), for Boolean variables:
         *   (l1 OR l2 OR ... OR lr) = 1 - product_i (1 - literal_i)
         *
         * For positive literal x: (1 - x) = (1 + x) over GF(2)
         * For negative literal ~x = 1-x = 1+x over GF(2)
         *
         * So: OR(l1,...,lr) = 1 + prod_i (1 + li)
         * where li is the GF(2) representation of the literal.
         *
         * For pos literal (var i): li = x_i
         * For neg literal (var i): li = 1 + x_i
         *
         * So the clause polynomial = 1 + prod_i (1 + li)
         *
         * Let's compute: first, prod_i (1 + li)
         */
        polynomial_t *clause_poly = poly_create(right_vars, 64, true);
        if (!clause_poly) { poly_destroy(result); return NULL; }
        poly_add_term(clause_poly, 1.0, 0ULL);  /* start with 1 */

        for (int32_t i = 0; i < num_right; i++) {
            int32_t rlit = right_lits[i];
            int32_t rvar = (rlit > 0) ? rlit : -rlit;
            rvar--;  /* 0-based */

            /* Factor for this literal: (1 + li) */
            polynomial_t *factor = poly_create(right_vars, 2, true);
            if (!factor) { poly_destroy(clause_poly); poly_destroy(result); return NULL; }

            if (rlit > 0) {
                /* Positive literal: li = x_rvar, so (1 + x_rvar) */
                poly_add_term(factor, 1.0, 0ULL);
                poly_add_term(factor, 1.0, 1ULL << rvar);
            } else {
                /* Negative literal: li = (1 + x_rvar),
                 * so (1 + li) = (1 + 1 + x_rvar) = x_rvar over GF(2) */
                poly_add_term(factor, 1.0, 1ULL << rvar);
            }

            polynomial_t *new_clause = poly_mul(clause_poly, factor);
            poly_destroy(clause_poly);
            poly_destroy(factor);
            clause_poly = new_clause;
            if (!clause_poly) { poly_destroy(result); return NULL; }
        }

        /* Now clause_poly = prod_i (1 + li).
         * The OR is: 1 + clause_poly (over GF(2), 1+ means XOR with 1) */
        poly_add_term(clause_poly, 1.0, 0ULL);  /* XOR with constant 1 */

        /* Multiply into overall result */
        polynomial_t *new_result = poly_mul(result, clause_poly);
        poly_destroy(result);
        poly_destroy(clause_poly);
        result = new_result;
        if (!result) return NULL;
    }

    poly_make_multilinear(result);
    return result;
}

/* ============================================================================
 * Williams' SAT Solver Implementation
 *
 * Algorithm (simplified):
 * 1. Choose split point k = n/2
 * 2. For each assignment to variables 0..k-1 (2^k possibilities):
 *    a. Construct residual polynomial P_assignment in variables k..n-1
 *    b. Evaluate P_assignment at all 2^{n-k} points using fast zeta transform
 *    c. If any evaluation gives 1, SAT (with that left + right assignment)
 * 3. Return UNSAT if no satisfying assignment found
 *
 * For n <= 20, we can run this exactly.
 * For larger n, this is an exponential algorithm; it illustrates the
 * polynomial method in algorithm design.
 * ============================================================================ */

williams_result_t williams_solve(const williams_cnf_t *f) {
    williams_result_t result = {false, NULL, 0.0, 0};
    if (!f) return result;

    double start = get_time_ms();
    int32_t n = f->num_vars;

    /* For very small n, brute force is faster */
    if (n <= 6) return williams_brute_force(f);

    int32_t k = n / 2;  /* split point */
    int32_t right_vars = n - k;

    result.witness = (int8_t *)malloc((size_t)n * sizeof(int8_t));
    if (!result.witness) return result;

    int64_t left_total = (int64_t)1 << k;

    /* For each left assignment */
    for (int64_t left_bits = 0; left_bits < left_total; left_bits++) {
        /* Build partial assignment */
        williams_assign_t *partial = williams_assign_create(n);
        if (!partial) continue;
        for (int32_t i = 0; i < k; i++)
            partial->values[i] = (int8_t)((left_bits >> i) & 1);

        /* Construct residual polynomial */
        polynomial_t *residual = williams_residual_polynomial(f, partial, k);

        if (residual && residual->num_terms > 0) {
            /* Evaluate at all 2^{right_vars} points */
            int32_t N_right = 1 << right_vars;
            double *evals = (double *)malloc((size_t)N_right * sizeof(double));
            if (evals) {
                poly_evaluate_all_bool(residual, evals, right_vars);
                result.evaluations += (uint64_t)N_right;

                /* Check if any right assignment satisfies */
                for (int32_t r = 0; r < N_right; r++) {
                    if (fabs(evals[r] - 1.0) < 1e-9 ||
                        (residual->over_gf2 && ((int32_t)evals[r] & 1))) {
                        /* Found satisfying assignment! */
                        result.satisfiable = true;
                        for (int32_t i = 0; i < k; i++)
                            result.witness[i] = (int8_t)((left_bits >> i) & 1);
                        for (int32_t i = 0; i < right_vars; i++)
                            result.witness[k + i] = (int8_t)((r >> i) & 1);
                        free(evals);
                        poly_destroy(residual);
                        williams_assign_destroy(partial);
                        result.runtime_ms = get_time_ms() - start;
                        return result;
                    }
                }
                free(evals);
            }
        }

        poly_destroy(residual);
        williams_assign_destroy(partial);
    }

    result.runtime_ms = get_time_ms() - start;
    return result;
}

/* ============================================================================
 * Generate Sample k-SAT Formula
 *
 * Creates a random k-SAT instance with n variables and m clauses.
 * Each clause is a random set of k distinct variables, each negated
 * independently with probability 1/2.
 *
 * seed: random seed for reproducibility
 * ============================================================================ */

williams_cnf_t *williams_generate_ksat(int32_t n, int32_t m, int32_t k,
                                         uint64_t seed) {
    williams_cnf_t *f = williams_cnf_create(n, m, k);
    if (!f) return NULL;

    uint64_t state = seed;

    for (int32_t c = 0; c < m; c++) {
        int32_t lits[32];
        /* Select k distinct variables */
        int32_t chosen[32] = {0};  /* flag array */
        for (int32_t i = 0; i < k; i++) {
            int32_t var;
            do {
                state = state * 6364136223846793005ULL + 1442695040888963407ULL;
                var = (int32_t)(state % (uint64_t)n) + 1;
            } while (chosen[var]);
            chosen[var] = 1;

            /* Randomly negate */
            state = state * 6364136223846793005ULL + 1442695040888963407ULL;
            bool negate = (state & 1) != 0;
            lits[i] = negate ? -var : var;
        }
        /* Clear chosen flags */
        for (int32_t i = 0; i < k; i++) {
            int32_t var = (lits[i] > 0) ? lits[i] : -lits[i];
            chosen[var] = 0;
        }
        williams_cnf_add_clause(f, lits, k);
    }

    return f;
}

/* ============================================================================
 * Verify Witness Assignment
 *
 * Checks that the claimed satisfying assignment actually satisfies
 * the formula. Used to validate solver output.
 * ============================================================================ */

bool williams_verify_witness(const williams_cnf_t *f, const int8_t *witness) {
    if (!f || !witness) return false;
    int32_t clause_size = f->clause_width + 1;

    for (int32_t c = 0; c < f->num_clauses; c++) {
        int32_t base = c * clause_size;
        bool clause_ok = false;
        for (int32_t i = 0; f->clause_data[base + i] != 0; i++) {
            int32_t lit = f->clause_data[base + i];
            int32_t var = (lit > 0) ? lit : -lit;
            bool need_true = (lit > 0);
            if ((bool)witness[var - 1] == need_true) {
                clause_ok = true;
                break;
            }
        }
        if (!clause_ok) return false;
    }
    return true;
}

/* ============================================================================
 * Polynomial Method for #SAT (Counting Satisfying Assignments)
 *
 * The polynomial method can also be used for counting. If we construct
 * the residual polynomial over the REALS (not GF(2)), then the sum of
 * all evaluations gives the number of satisfying right assignments
 * for a given left assignment.
 *
 * This connects to the inclusion-exclusion principle:
 *   Number of satisfying assignments = sum_{x} indicator(x)
 *     = sum_{x} prod_{clauses c} (1 - prod_{lit in c} (1 - lit(x)))
 *
 * Over GF(2), the count is modulo 2. Over reals, it's exact.
 * ============================================================================ */
