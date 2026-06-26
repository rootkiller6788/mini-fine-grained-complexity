#ifndef POLY_METHOD_H
#define POLY_METHOD_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <math.h>

/* ============================================================================
 * mini-polynomial-method-consequences
 *
 * Polynomial method applications to fine-grained complexity.
 * Core definitions, data structures, and algorithms for using polynomial
 * representations to achieve faster-than-brute-force algorithms.
 *
 * Key references:
 *   - Razborov (1987): "Lower bounds on the size of bounded-depth circuits
 *     over a complete basis with logical addition"
 *   - Smolensky (1987): "Algebraic methods in the theory of lower bounds
 *     for Boolean circuit complexity"
 *   - Williams (2014): "Nonuniform ACC Circuit Lower Bounds"
 *   - Abboud, Williams, Yu (2015): "More Applications of the Polynomial
 *     Method to Algorithm Design"
 *
 * L1: Core definitions (polynomial, monomial, algebraic circuit)
 * L2: Core concepts (degree, modulus switching, polynomial representation)
 * L3: Mathematical structures (GF(p) arithmetic, multi-linear extensions)
 * L4: Fundamental laws (Razborov-Smolensky, Williams' algorithm correctness)
 * ============================================================================ */

/* ---- GF(2) Element ---- */
typedef uint8_t gf2_t;

/* ---- A monomial: product of variables, represented as a bitmask ---- */
typedef struct {
    uint64_t vars;       /* bitmask: bit i set iff variable i appears */
    int32_t  degree;     /* number of variables in the monomial */
} monomial_t;

/* ---- A term: coefficient * monomial ---- */
typedef struct {
    double     coeff;    /* coefficient (0.0 means term absent) */
    monomial_t mono;     /* the monomial */
} term_t;

/* ---- A multivariate polynomial over GF(2) or reals, sparse representation ---- */
typedef struct {
    term_t    *terms;     /* array of terms */
    int32_t    num_terms; /* current number of terms */
    int32_t    capacity;  /* allocated capacity */
    int32_t    num_vars;  /* number of variables */
    int32_t    degree;    /* total degree of polynomial */
    bool       over_gf2;  /* true if coefficients are in GF(2) */
} polynomial_t;

/* ---- A multi-linear polynomial ---- */
typedef struct {
    polynomial_t *poly;  /* underlying polynomial, guaranteed multi-linear */
} multilinear_poly_t;

/* ---- A degree-d polynomial approximation to a Boolean function ---- */
typedef struct {
    polynomial_t *approximator;
    int32_t       num_vars;
    int32_t       degree;
    double        error;
} poly_approximation_t;

/* ---- An algebraic circuit over a field F ---- */
typedef enum { AC_INPUT, AC_CONSTANT, AC_ADD, AC_MUL } ac_gate_type_t;

typedef struct {
    ac_gate_type_t type;
    int32_t        var_idx;
    double         constant;
    int32_t        left;
    int32_t        right;
} ac_gate_t;

typedef struct {
    ac_gate_t *gates;
    int32_t    num_gates;
    int32_t    capacity;
    int32_t    num_vars;
    int32_t    output_gate;
} algebraic_circuit_t;

/* ---- The Orthogonal Vectors (OV) problem instance ---- */
typedef struct {
    uint64_t *vectors;
    int32_t   num_vectors;
    int32_t   dimension;
} ov_instance_t;

/* ---- OV result ---- */
typedef struct {
    bool    orthogonal_pair_exists;
    int32_t a_idx;
    int32_t b_idx;
    int32_t dot_product;
} ov_result_t;

/* ---- CNF formula for Williams' algorithm ---- */
typedef struct {
    int32_t *clause_data;
    int32_t  num_clauses;
    int32_t  num_vars;
    int32_t  clause_width;
    int32_t  data_len;
} williams_cnf_t;

/* ---- Partial assignment for Williams' splitting ---- */
typedef struct {
    int8_t  *values;
    int32_t  num_vars;
} williams_assign_t;

/* ---- Result of Williams-style satisfiability check ---- */
typedef struct {
    bool      satisfiable;
    int8_t   *witness;
    double    runtime_ms;
    uint64_t  evaluations;
} williams_result_t;

/* ============================================================================
 * GF(p) Arithmetic and Polynomial Operations
 * ============================================================================ */

static inline gf2_t gf2_add(gf2_t a, gf2_t b)  { return a ^ b; }
static inline gf2_t gf2_mul(gf2_t a, gf2_t b)  { return a & b; }
static inline gf2_t gf2_sub(gf2_t a, gf2_t b)  { return a ^ b; }

typedef struct { int32_t prime; } gf_p_t;

int32_t gf_p_add(int32_t a, int32_t b, int32_t p);
int32_t gf_p_mul(int32_t a, int32_t b, int32_t p);
int32_t gf_p_pow(int32_t base, int32_t exp, int32_t p);
int32_t gf_p_neg(int32_t a, int32_t p);
int32_t gf_p_inv(int32_t a, int32_t p);

/* Polynomial lifecycle */
polynomial_t *poly_create(int32_t num_vars, int32_t capacity, bool over_gf2);
void          poly_destroy(polynomial_t *p);
polynomial_t *poly_clone(const polynomial_t *p);

/* Term manipulation */
void poly_add_term(polynomial_t *p, double coeff, uint64_t var_mask);
void poly_add_term_int(polynomial_t *p, int32_t coeff, uint64_t var_mask);

/* Polynomial evaluation */
double  poly_evaluate_bool(const polynomial_t *p, uint64_t x);
int32_t poly_evaluate_gf2(const polynomial_t *p, uint64_t x);

/* Polynomial arithmetic */
polynomial_t *poly_add(const polynomial_t *a, const polynomial_t *b);
polynomial_t *poly_mul(const polynomial_t *a, const polynomial_t *b);
polynomial_t *poly_scalar_mul(const polynomial_t *p, double c);

/* Multi-linear reduction */
void poly_make_multilinear(polynomial_t *p);

/* Degree computation */
int32_t poly_degree(const polynomial_t *p);

/* Debug printing */
void poly_print(const polynomial_t *p);

/* ============================================================================
 * Polynomial Representations of Boolean Functions
 * ============================================================================ */

/* Truth table to GF(2) polynomial via Mobius transform (Zhegalkin/ANF).
 * Theorem: Every Boolean f has a unique multi-linear polynomial over GF(2).
 * Time: O(n * 2^n). */
polynomial_t *truth_table_to_poly_gf2(const uint8_t *truth_table, int32_t n);

/* Truth table to Fourier polynomial over reals.
 * Time: O(n * 2^n). */
polynomial_t *truth_table_to_poly_fourier(const uint8_t *truth_table, int32_t n);

/* Fourier coefficient f_hat(S) for subset S */
double fourier_coefficient(const uint8_t *truth_table, int32_t n, uint64_t subset);

/* ============================================================================
 * Fast Transforms and Evaluation
 * ============================================================================ */

/* Fast multi-point evaluation at all 2^n Boolean points.
 * Uses the fast subset sum / zeta transform. Time: O(n * 2^n). */
void poly_evaluate_all_bool(const polynomial_t *p, double *output, int32_t n);

/* Fast Walsh-Hadamard Transform. Time: O(n * 2^n). */
void walsh_hadamard_transform(double *data, int32_t n);

/* Polynomial AND: pointwise product of two Boolean functions. */
polynomial_t *poly_and(const polynomial_t *f, const polynomial_t *g);

/* Polynomial XOR: pointwise XOR over GF(2). */
polynomial_t *poly_xor(const polynomial_t *f, const polynomial_t *g);

/* ============================================================================
 * OV Problem Solvers
 * ============================================================================ */

ov_instance_t *ov_create(int32_t n, int32_t d);
void           ov_destroy(ov_instance_t *ov);
void           ov_set_vector(ov_instance_t *ov, int32_t idx, uint64_t vec);
uint64_t       ov_get_vector(const ov_instance_t *ov, int32_t idx);

/* Brute-force OV: O(n^2 * d) */
ov_result_t ov_solve_brute_force(const ov_instance_t *A, const ov_instance_t *B);

/* Polynomial-method OV: O(n^{2-epsilon}) for d = O(log n) */
ov_result_t ov_solve_polynomial(const ov_instance_t *A, const ov_instance_t *B);

/* Word-packed OV: O(n^2) with fast popcount, d <= 64 */
ov_result_t ov_solve_packed(const ov_instance_t *A, const ov_instance_t *B);

/* Count orthogonal pairs */
int64_t ov_count_orthogonal_pairs(const ov_instance_t *A, const ov_instance_t *B);

/* ============================================================================
 * Williams Algorithm for CNF-SAT
 * ============================================================================ */

williams_cnf_t   *williams_cnf_create(int32_t num_vars, int32_t num_clauses, int32_t k);
void              williams_cnf_destroy(williams_cnf_t *f);
void              williams_cnf_add_clause(williams_cnf_t *f, const int32_t *lits, int32_t len);
williams_assign_t *williams_assign_create(int32_t num_vars);
void               williams_assign_destroy(williams_assign_t *a);
bool               williams_cnf_evaluate(const williams_cnf_t *f, const williams_assign_t *a);

/* Williams-style SAT solver using polynomial method.
 * Splits n variables -> O*(2^{n/2}) time. */
williams_result_t williams_solve(const williams_cnf_t *f);

/* Residual polynomial for partial assignment */
polynomial_t *williams_residual_polynomial(const williams_cnf_t *f,
                                             const williams_assign_t *partial,
                                             int32_t split_point);

/* Brute-force SAT baseline */
williams_result_t williams_brute_force(const williams_cnf_t *f);

/* ============================================================================
 * Additional Polynomial Utilities
 * ============================================================================ */

/* Inverse Walsh-Hadamard transform */
void inverse_walsh_hadamard_transform(double *data, int32_t n);

/* Polynomial degree spectrum: mass at each degree level */
double *poly_degree_spectrum(const polynomial_t *p, int32_t n);

/* Polynomial sparsity: number of non-zero terms */
int32_t poly_sparsity(const polynomial_t *p);

/* Check if polynomial is multi-linear */
bool poly_is_multilinear(const polynomial_t *p);

/* ============================================================================
 * OV Additional Functions
 * ============================================================================ */

/* Generate random OV instance with given seed */
ov_instance_t *ov_generate_random(int32_t n, int32_t d, uint64_t seed);

/* Verify a claimed orthogonal pair */
bool ov_verify_orthogonal(const ov_instance_t *A, const ov_instance_t *B,
                           int32_t a_idx, int32_t b_idx);

/* Polynomial method OV with batch evaluation */
ov_result_t ov_solve_polynomial_batch(const ov_instance_t *A,
                                        const ov_instance_t *B);

/* Polynomial method OV using zeta transform */
ov_result_t ov_solve_polynomial_zeta(const ov_instance_t *A,
                                       const ov_instance_t *B);

/* ============================================================================
 * Williams Algorithm Additional Functions
 * ============================================================================ */

/* Generate random k-SAT formula */
williams_cnf_t *williams_generate_ksat(int32_t n, int32_t m, int32_t k,
                                         uint64_t seed);

/* Verify a witness assignment */
bool williams_verify_witness(const williams_cnf_t *f, const int8_t *witness);

/* ============================================================================
 * Algebraic Circuit Operations (L8)
 * ============================================================================ */

algebraic_circuit_t *ac_create(int32_t capacity, int32_t num_vars);
void ac_destroy(algebraic_circuit_t *circ);
int32_t ac_add_input(algebraic_circuit_t *circ, int32_t var_idx);
int32_t ac_add_constant(algebraic_circuit_t *circ, double value);
int32_t ac_add_add(algebraic_circuit_t *circ, int32_t left, int32_t right);
int32_t ac_add_mul(algebraic_circuit_t *circ, int32_t left, int32_t right);
double ac_evaluate(const algebraic_circuit_t *circ, const double *input_values);
int32_t ac_compute_degree(const algebraic_circuit_t *circ);
polynomial_t *ac_to_polynomial(const algebraic_circuit_t *circ);
algebraic_circuit_t *ac_build_mod_p(int32_t n, int32_t p);
algebraic_circuit_t *ac_build_parity(int32_t n);

/* Approximate AND of k variables over GF(p) */
polynomial_t *ac_approximate_and(int32_t k, int32_t p, int32_t degree);

/* ============================================================================
 * Boolean Function Analysis (L4, L8)
 * ============================================================================ */

/* Build truth table for MOD_m function */
uint8_t *build_mod_m_truth_table(int32_t n, int32_t m);

/* Build truth table for MAJORITY function */
uint8_t *build_majority_truth_table(int32_t n);

/* Compute GF(p)-degree of a Boolean function from its truth table */
int32_t gf_degree_from_truth_table(const uint8_t *truth_table, int32_t n, int32_t p);

/* Compute degree of MOD_m over GF(p) */
int32_t mod_m_degree_over_p(int32_t n, int32_t m, int32_t p);

/* Compute degree of MAJORITY over GF(p) */
int32_t majority_degree_over_p(int32_t n, int32_t p);

/* Compute decision tree degree lower bound */
int32_t decision_tree_degree_lower_bound(const uint8_t *truth_table, int32_t n);

/* Check if function has low-degree concentration (epsilon up to degree d) */
bool function_has_low_degree_concentration(const uint8_t *truth_table,
                                             int32_t n, int32_t d,
                                             double epsilon);

/* Compute variable influences */
double *compute_influences(const uint8_t *truth_table, int32_t n);

/* Boolean extension degree */
int32_t boolean_extension_degree(const uint8_t *truth_table, int32_t n, int32_t p);

/* ============================================================================
 * Razborov-Smolensky Framework (L4, L8)
 * ============================================================================ */

/* Parameters for Razborov approximation analysis */
typedef struct {
    int32_t depth;
    int32_t size;
    int32_t max_fan_in;
    int32_t p;
} razborov_approx_params_t;

/* Compute polynomial degree bound for given circuit parameters */
int32_t razborov_degree_bound(const razborov_approx_params_t *params,
                                double epsilon);

/* Check if a function's degree exceeds the Razborov bound */
bool razborov_lower_bound_check(int32_t function_degree,
                                  const razborov_approx_params_t *params,
                                  double epsilon);

/* Result of Smolensky analysis */
typedef struct {
    int32_t n;
    int32_t mod_m;
    int32_t field_p;
    int32_t exact_degree;
    bool requires_high_degree;
} smolensky_result_t;

/* Analyze MOD_m function over GF(p) */
smolensky_result_t smolensky_analyze(int32_t n, int32_t mod_m, int32_t field_p);

/* Natural property check */
typedef struct {
    bool is_large;
    bool is_useful;
    double largeness;
    int32_t degree_threshold;
} natural_property_check_t;

natural_property_check_t analyze_truth_table_property(
    const uint8_t *truth_table, int32_t n);

/* Circuit complexity estimate */
typedef struct {
    int32_t n;
    int32_t truth_table_size;
    int32_t exact_degree;
    int32_t sparsity;
    double complexity_estimate;
} circuit_complexity_estimate_t;

circuit_complexity_estimate_t estimate_circuit_complexity(
    const uint8_t *truth_table, int32_t n);

/* Degree comparison for MOD functions */
typedef struct {
    int32_t n;
    int32_t mod_m;
    int32_t field_p;
    int32_t degree;
    double ratio;
} degree_comparison_t;

degree_comparison_t compare_mod_degrees(int32_t n);

#endif /* POLY_METHOD_H */
