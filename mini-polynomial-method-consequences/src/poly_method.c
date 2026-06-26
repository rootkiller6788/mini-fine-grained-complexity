#include "poly_method.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ============================================================================
 * poly_method.c - Core Polynomial Operations
 *
 * L3: Mathematical structures - GF(p) arithmetic, polynomial rings
 * L4: Fundamental laws - polynomial representation theorems
 * L5: Algorithms - fast transforms (FWHT, zeta), evaluation, arithmetic
 * ============================================================================ */

/* ---- GF(p) Modular Arithmetic ---- */

int32_t gf_p_add(int32_t a, int32_t b, int32_t p) {
    int32_t r = (a + b) % p;
    return (r < 0) ? r + p : r;
}

int32_t gf_p_mul(int32_t a, int32_t b, int32_t p) {
    int32_t r = (int32_t)(((int64_t)a * (int64_t)b) % p);
    return (r < 0) ? r + p : r;
}

int32_t gf_p_pow(int32_t base, int32_t exp, int32_t p) {
    int32_t result = 1;
    int32_t b = ((base % p) + p) % p;
    int32_t e = exp;
    while (e > 0) {
        if (e & 1) result = gf_p_mul(result, b, p);
        b = gf_p_mul(b, b, p);
        e >>= 1;
    }
    return result;
}

int32_t gf_p_neg(int32_t a, int32_t p) {
    if (a == 0) return 0;
    int32_t r = a % p;
    if (r < 0) r += p;
    return (r == 0) ? 0 : p - r;
}

/* Extended Euclidean algorithm for modular inverse */
int32_t gf_p_inv(int32_t a, int32_t p) {
    int32_t t = 0, newt = 1;
    int32_t r = p, newr = a % p;
    if (newr < 0) newr += p;
    if (newr == 0) return 0;
    while (newr != 0) {
        int32_t q = r / newr;
        int32_t tmp = newt;
        newt = t - q * newt;
        t = tmp;
        tmp = newr;
        newr = r - q * newr;
        r = tmp;
    }
    if (r > 1) return 0;
    if (t < 0) t += p;
    return t;
}

/* ---- Polynomial Lifecycle ---- */

polynomial_t *poly_create(int32_t num_vars, int32_t capacity, bool over_gf2) {
    polynomial_t *p = (polynomial_t *)malloc(sizeof(polynomial_t));
    if (!p) return NULL;
    if (capacity <= 0) capacity = 64;
    p->terms = (term_t *)calloc((size_t)capacity, sizeof(term_t));
    if (!p->terms) { free(p); return NULL; }
    p->num_terms = 0;
    p->capacity = capacity;
    p->num_vars = num_vars;
    p->degree = 0;
    p->over_gf2 = over_gf2;
    return p;
}

void poly_destroy(polynomial_t *p) {
    if (p) { free(p->terms); free(p); }
}

polynomial_t *poly_clone(const polynomial_t *p) {
    if (!p) return NULL;
    polynomial_t *c = poly_create(p->num_vars, p->capacity, p->over_gf2);
    if (!c) return NULL;
    c->num_terms = p->num_terms;
    c->degree = p->degree;
    memcpy(c->terms, p->terms, (size_t)p->num_terms * sizeof(term_t));
    return c;
}

/* ---- Term Manipulation ---- */

static int32_t popcount64(uint64_t x) {
    x = x - ((x >> 1) & 0x5555555555555555ULL);
    x = (x & 0x3333333333333333ULL) + ((x >> 2) & 0x3333333333333333ULL);
    x = (x + (x >> 4)) & 0x0F0F0F0F0F0F0F0FULL;
    return (int32_t)((x * 0x0101010101010101ULL) >> 56);
}

void poly_add_term(polynomial_t *p, double coeff, uint64_t var_mask) {
    if (!p || coeff == 0.0) return;
    if (p->over_gf2) {
        int32_t ci = (int32_t)coeff;
        coeff = (double)(ci & 1);
        if (coeff == 0.0) return;
    }
    /* Check existing term with same variable mask */
    for (int32_t i = 0; i < p->num_terms; i++) {
        if (p->terms[i].mono.vars == var_mask && p->terms[i].coeff != 0.0) {
            if (p->over_gf2) {
                int32_t nc = ((int32_t)p->terms[i].coeff) ^ ((int32_t)coeff);
                p->terms[i].coeff = (double)nc;
                if (nc == 0) {
                    p->terms[i] = p->terms[p->num_terms - 1];
                    p->num_terms--;
                    int32_t md = 0;
                    for (int32_t j = 0; j < p->num_terms; j++)
                        if (p->terms[j].coeff != 0.0 && p->terms[j].mono.degree > md)
                            md = p->terms[j].mono.degree;
                    p->degree = md;
                }
            } else {
                p->terms[i].coeff += coeff;
                if (fabs(p->terms[i].coeff) < 1e-12) {
                    p->terms[i] = p->terms[p->num_terms - 1];
                    p->num_terms--;
                    int32_t md = 0;
                    for (int32_t j = 0; j < p->num_terms; j++)
                        if (fabs(p->terms[j].coeff) > 1e-12 && p->terms[j].mono.degree > md)
                            md = p->terms[j].mono.degree;
                    p->degree = md;
                }
            }
            return;
        }
    }
    /* Grow if needed */
    if (p->num_terms >= p->capacity) {
        int32_t nc = p->capacity * 2;
        term_t *nt = (term_t *)realloc(p->terms, (size_t)nc * sizeof(term_t));
        if (!nt) return;
        p->terms = nt;
        p->capacity = nc;
    }
    int32_t deg = popcount64(var_mask);
    p->terms[p->num_terms].coeff = coeff;
    p->terms[p->num_terms].mono.vars = var_mask;
    p->terms[p->num_terms].mono.degree = deg;
    p->num_terms++;
    if (deg > p->degree) p->degree = deg;
}

void poly_add_term_int(polynomial_t *p, int32_t coeff, uint64_t var_mask) {
    poly_add_term(p, (double)coeff, var_mask);
}

/* ---- Polynomial Evaluation ---- */

double poly_evaluate_bool(const polynomial_t *p, uint64_t x) {
    if (!p) return 0.0;
    double result = 0.0;
    for (int32_t i = 0; i < p->num_terms; i++) {
        if (fabs(p->terms[i].coeff) < 1e-12) continue;
        if ((p->terms[i].mono.vars & x) == p->terms[i].mono.vars)
            result += p->terms[i].coeff;
    }
    if (p->over_gf2) return (double)(((int32_t)result) & 1);
    return result;
}

int32_t poly_evaluate_gf2(const polynomial_t *p, uint64_t x) {
    if (!p) return 0;
    int32_t result = 0;
    for (int32_t i = 0; i < p->num_terms; i++) {
        if (fabs(p->terms[i].coeff) < 1e-12) continue;
        if ((p->terms[i].mono.vars & x) == p->terms[i].mono.vars)
            result ^= ((int32_t)p->terms[i].coeff) & 1;
    }
    return result;
}

/* ---- Polynomial Arithmetic ---- */

polynomial_t *poly_add(const polynomial_t *a, const polynomial_t *b) {
    if (!a || !b) return NULL;
    int32_t mv = (a->num_vars > b->num_vars) ? a->num_vars : b->num_vars;
    polynomial_t *r = poly_create(mv, a->num_terms + b->num_terms,
                                   a->over_gf2 && b->over_gf2);
    if (!r) return NULL;
    for (int32_t i = 0; i < a->num_terms; i++)
        poly_add_term(r, a->terms[i].coeff, a->terms[i].mono.vars);
    for (int32_t i = 0; i < b->num_terms; i++)
        poly_add_term(r, b->terms[i].coeff, b->terms[i].mono.vars);
    return r;
}

polynomial_t *poly_mul(const polynomial_t *a, const polynomial_t *b) {
    if (!a || !b) return NULL;
    int32_t mv = (a->num_vars > b->num_vars) ? a->num_vars : b->num_vars;
    bool gf2 = a->over_gf2 && b->over_gf2;
    polynomial_t *r = poly_create(mv, a->num_terms * b->num_terms, gf2);
    if (!r) return NULL;
    for (int32_t i = 0; i < a->num_terms; i++) {
        if (fabs(a->terms[i].coeff) < 1e-12) continue;
        for (int32_t j = 0; j < b->num_terms; j++) {
            if (fabs(b->terms[j].coeff) < 1e-12) continue;
            double nc;
            if (gf2)
                nc = (double)(((int32_t)a->terms[i].coeff & 1) &
                              ((int32_t)b->terms[j].coeff & 1));
            else
                nc = a->terms[i].coeff * b->terms[j].coeff;
            uint64_t nv = a->terms[i].mono.vars | b->terms[j].mono.vars;
            poly_add_term(r, nc, nv);
        }
    }
    return r;
}

polynomial_t *poly_scalar_mul(const polynomial_t *p, double c) {
    if (!p || fabs(c) < 1e-12) return NULL;
    polynomial_t *r = poly_clone(p);
    if (!r) return NULL;
    if (p->over_gf2) {
        if (((int32_t)c & 1) == 0) { r->num_terms = 0; r->degree = 0; }
        return r;
    }
    for (int32_t i = 0; i < r->num_terms; i++)
        r->terms[i].coeff *= c;
    return r;
}

/* ---- Multi-linear Reduction ---- */

void poly_make_multilinear(polynomial_t *p) {
    if (!p) return;
    for (int32_t i = 0; i < p->num_terms; i++) {
        int32_t d = popcount64(p->terms[i].mono.vars);
        p->terms[i].mono.degree = (d > p->num_vars) ? p->num_vars : d;
    }
    int32_t md = 0;
    for (int32_t i = 0; i < p->num_terms; i++)
        if (fabs(p->terms[i].coeff) > 1e-12 && p->terms[i].mono.degree > md)
            md = p->terms[i].mono.degree;
    p->degree = md;
}

int32_t poly_degree(const polynomial_t *p) {
    return p ? p->degree : -1;
}

/* ---- Debug Print ---- */

void poly_print(const polynomial_t *p) {
    if (!p) { printf("NULL\n"); return; }
    printf("Poly(vars=%d, terms=%d, deg=%d, gf2=%d):\n  ",
           p->num_vars, p->num_terms, p->degree, p->over_gf2);
    if (p->num_terms == 0) { printf("0\n"); return; }
    int printed = 0;
    for (int32_t i = 0; i < p->num_terms; i++) {
        if (p->over_gf2 && ((int32_t)p->terms[i].coeff & 1) == 0) continue;
        if (!p->over_gf2 && fabs(p->terms[i].coeff) < 1e-12) continue;
        if (printed > 0) printf(" + ");
        if (p->over_gf2)
            printf("%d", (int)p->terms[i].coeff & 1);
        else
            printf("%.3f", p->terms[i].coeff);
        uint64_t v = p->terms[i].mono.vars;
        int32_t bit = 0;
        int started = 0;
        while (v) {
            if (v & 1) { printf("*x%d", bit); started = 1; }
            v >>= 1;
            bit++;
        }
        if (!started && printed == 0 && p->terms[i].mono.vars == 0) {
            /* constant term: coefficient already printed */
        }
        printed++;
    }
    printf("\n");
}

/* ============================================================================
 * Fast Walsh-Hadamard Transform (FWHT)
 *
 * Computes the Fourier transform of f: {0,1}^n -> R over (Z_2)^n.
 * For f: {0,1}^n -> R, the transform gives:
 *   F_hat[S] = sum_{x in {0,1}^n} f(x) * (-1)^{<S,x>}
 *
 * where <S,x> = sum_{i in S} x_i (mod 2) is the inner product.
 * This is the Fourier transform on the Boolean cube.
 *
 * Algorithm: In-place iterative butterfly, O(n * 2^n) time.
 * For n variables, array size must be 2^n.
 *
 * Reference: O'Donnell, "Analysis of Boolean Functions" (2014), Chapter 1.
 *
 * Applications in complexity:
 * - Computing the Fourier spectrum of Boolean functions
 * - Learning parities (Goldreich-Levin theorem)
 * - Proving circuit lower bounds via Fourier concentration
 * - Testing whether a function is a junta
 * ============================================================================ */

void walsh_hadamard_transform(double *data, int32_t n) {
    int32_t N = 1 << n;
    if (N <= 0 || !data) return;
    for (int32_t step = 1; step < N; step <<= 1) {
        for (int32_t i = 0; i < N; i += (step << 1)) {
            for (int32_t j = i; j < i + step; j++) {
                double u = data[j];
                double v = data[j + step];
                data[j] = u + v;
                data[j + step] = u - v;
            }
        }
    }
}

/* ============================================================================
 * Inverse Walsh-Hadamard Transform
 *
 * Since FWHT is self-inverse up to scaling:
 *   IWHT(data) = FWHT(data) / 2^n
 *
 * Time: O(n * 2^n).
 * ============================================================================ */

void inverse_walsh_hadamard_transform(double *data, int32_t n) {
    walsh_hadamard_transform(data, n);
    int32_t N = 1 << n;
    double scale = 1.0 / (double)N;
    for (int32_t i = 0; i < N; i++) data[i] *= scale;
}

/* ============================================================================
 * Zeta Transform (Subset Sum / Superset Sum) and Mobius Transform
 *
 * Zeta transform (subset): (Z f)[T] = sum_{S subseteq T} f(S)
 * Zeta transform (superset): (Z' f)[T] = sum_{S supseteq T} f(S)
 * Mobius transform (subset): (M g)[S] = sum_{T subseteq S} (-1)^{|S|-|T|} g(T)
 *
 * These are mutual inverses. The zeta transform on the coefficient vector
 * of a multi-linear polynomial directly computes the polynomial's value
 * at all Boolean points. This is the core insight of the polynomial method:
 * evaluating a polynomial at all 2^n points costs O(n * 2^n), which is
 * linear in the output size.
 *
 * Time: O(n * 2^n).
 *
 * Reference: Bjorklund, Husfeldt, Kaski, Koivisto (2007),
 *   "Fourier meets Mobius: fast subset convolution".
 * ============================================================================ */

/* Generic zeta-like transform: sign=+1 for zeta, sign=-1 for mobius */
static void zeta_transform_generic(double *f, int32_t n, double sign) {
    int32_t N = 1 << n;
    for (int32_t i = 0; i < n; i++) {
        for (int32_t mask = 0; mask < N; mask++) {
            if (mask & (1 << i)) {
                f[mask] += sign * f[mask ^ (1 << i)];
            }
        }
    }
}

static void zeta_transform(double *f, int32_t n) {
    zeta_transform_generic(f, n, 1.0);
}

static void mobius_transform(double *f, int32_t n) {
    zeta_transform_generic(f, n, -1.0);
}

/* Superset zeta transform: g[T] = sum_{S supseteq T} f[S] */
static void zeta_superset_transform(double *f, int32_t n) {
    int32_t N = 1 << n;
    for (int32_t i = 0; i < n; i++) {
        for (int32_t mask = 0; mask < N; mask++) {
            if (!(mask & (1 << i))) {
                f[mask] += f[mask | (1 << i)];
            }
        }
    }
}

/* Superset Mobius transform (inverse of superset zeta) */
static void mobius_superset_transform(double *f, int32_t n) {
    int32_t N = 1 << n;
    for (int32_t i = 0; i < n; i++) {
        for (int32_t mask = 0; mask < N; mask++) {
            if (!(mask & (1 << i))) {
                f[mask] -= f[mask | (1 << i)];
            }
        }
    }
}

/* ============================================================================
 * Fast Multi-Point Boolean Evaluation
 *
 * For a multi-linear polynomial p(x) = sum_{S subseteq [n]} c_S * prod_{i in S} x_i,
 * evaluating at Boolean point x gives:
 *   p(x) = sum_{S subseteq x} c_S
 *
 * This is exactly the zeta transform of the coefficient vector:
 *   p_values = ZetaTransform(coefficients)
 *
 * So we can evaluate p at ALL 2^n Boolean points in O(n * 2^n) time
 * by computing the zeta transform of c.
 *
 * This is a KEY subroutine in the polynomial method for algorithm design:
 * - Williams' SAT algorithm evaluates 2^{n/2} polynomials at 2^{n/2} points
 *   in O*(2^{n/2}) time using this technique.
 * - The OV polynomial method evaluates structured polynomials at many points.
 * - Circuit analysis: evaluate approximating polynomials at all inputs.
 *
 * Time: O(n * 2^n). Output: 2^n values.
 * ============================================================================ */

void poly_evaluate_all_bool(const polynomial_t *p, double *output, int32_t n) {
    int32_t N = 1 << n;
    if (!p || !output || N <= 0) return;

    /* Initialize with coefficients c[S] in position S */
    memset(output, 0, (size_t)N * sizeof(double));
    for (int32_t i = 0; i < p->num_terms; i++) {
        if (fabs(p->terms[i].coeff) < 1e-12) continue;
        uint64_t vars = p->terms[i].mono.vars;
        if (vars < (uint64_t)N) {
            if (p->over_gf2) {
                output[vars] = (double)((int32_t)output[vars] ^
                                ((int32_t)p->terms[i].coeff & 1));
            } else {
                output[vars] += p->terms[i].coeff;
            }
        }
    }

    /* Zeta transform: output[x] = sum_{S subseteq x} c[S] = p(x) */
    zeta_transform(output, n);
}

/* ============================================================================
 * Truth Table to GF(2) Polynomial (Algebraic Normal Form / Zhegalkin Polynomial)
 *
 * Theorem (Zhegalkin 1927, Stone 1936): Every Boolean function
 * f: {0,1}^n -> {0,1} has a UNIQUE representation as a multi-linear
 * polynomial over GF(2):
 *
 *   f(x_1,...,x_n) = XOR_{S subseteq [n]} a_S * AND_{i in S} x_i
 *
 * The coefficients a_S are computed by the GF(2) Mobius transform:
 *   a_S = XOR_{T subseteq S} f(T)
 *
 * Proof sketch:
 *   Define g(S) = f(char_vector(S)). Then a_S = XOR_{T subseteq S} g(T).
 *   The inverse is g(T) = XOR_{S subseteq T} a_S, which is exactly
 *   the evaluation formula (since x_i = 1 iff i is in the assignment).
 *
 * This isomorphism between Boolean functions and GF(2) polynomials
 * is fundamental to:
 * - Circuit lower bounds (Razborov 1987, Smolensky 1987)
 * - Algebraic attacks on cryptography
 * - Property testing and learning theory
 * - The polynomial method in algorithm design
 *
 * Time: O(n * 2^n). Space: O(2^n).
 * ============================================================================ */

polynomial_t *truth_table_to_poly_gf2(const uint8_t *truth_table, int32_t n) {
    int32_t N = 1 << n;
    if (!truth_table || N <= 0 || n > 24) return NULL;

    double *coeffs = (double *)calloc((size_t)N, sizeof(double));
    if (!coeffs) return NULL;

    for (int32_t i = 0; i < N; i++)
        coeffs[i] = (double)(truth_table[i] & 1);

    /* GF(2) Mobius transform: a_S = XOR_{T subseteq S} f(T) */
    for (int32_t i = 0; i < n; i++) {
        for (int32_t mask = 0; mask < N; mask++) {
            if (mask & (1 << i)) {
                int32_t a = (int32_t)coeffs[mask];
                int32_t b = (int32_t)coeffs[mask ^ (1 << i)];
                coeffs[mask] = (double)(a ^ b);
            }
        }
    }

    polynomial_t *result = poly_create(n, N, true);
    for (int32_t S = 0; S < N; S++) {
        if ((int32_t)coeffs[S] & 1)
            poly_add_term(result, 1.0, (uint64_t)S);
    }

    free(coeffs);
    return result;
}

/* ============================================================================
 * Truth Table to Fourier Polynomial over Reals
 *
 * The Fourier expansion represents f: {0,1}^n -> R as:
 *   f(x) = sum_{S subseteq [n]} f_hat(S) * chi_S(x)
 *
 * where chi_S(x) = (-1)^{sum_{i in S} x_i} = prod_{i in S} (1 - 2*x_i)
 *
 * Properties of the Fourier basis:
 * - Orthonormal: <chi_S, chi_T> = delta_{S,T}
 * - Parseval: sum_S f_hat(S)^2 = E_x[f(x)^2]
 * - Convolution: (f*g)_hat(S) = sum_{T} f_hat(T) * g_hat(S Delta T)
 *
 * The Fourier coefficients are:
 *   f_hat(S) = E_x[f(x) * chi_S(x)] = (1/2^n) * sum_x f(x) * chi_S(x)
 *
 * The Fourier expansion is a powerful tool in:
 * - Analysis of Boolean functions (O'Donnell 2014)
 * - Hardness of approximation (Hastad's 3-bit PCP)
 * - Learning theory (low-degree algorithm, KM93)
 * - Quantum computing (phase polynomials)
 *
 * Time: O(n * 2^n) via FWHT.
 * ============================================================================ */

polynomial_t *truth_table_to_poly_fourier(const uint8_t *truth_table, int32_t n) {
    int32_t N = 1 << n;
    if (!truth_table || N <= 0 || n > 24) return NULL;

    double *f_vals = (double *)calloc((size_t)N, sizeof(double));
    if (!f_vals) return NULL;
    for (int32_t i = 0; i < N; i++)
        f_vals[i] = (double)(truth_table[i] & 1);

    /* FWHT: result[S] = sum_x f(x) * chi_S(x) */
    walsh_hadamard_transform(f_vals, n);

    /* Normalize: f_hat(S) = result[S] / 2^n */
    double norm = 1.0 / (double)N;
    for (int32_t i = 0; i < N; i++)
        f_vals[i] *= norm;

    polynomial_t *result = poly_create(n, N, false);
    for (int32_t S = 0; S < N; S++) {
        if (fabs(f_vals[S]) > 1e-15)
            poly_add_term(result, f_vals[S], (uint64_t)S);
    }

    free(f_vals);
    return result;
}

/* ============================================================================
 * Individual Fourier Coefficient
 *
 * f_hat(S) = (1/2^n) * sum_{x in {0,1}^n} f(x) * (-1)^{<S,x>}
 *
 * Computes a single coefficient in O(2^n) without full transform.
 * ============================================================================ */

double fourier_coefficient(const uint8_t *truth_table, int32_t n, uint64_t subset) {
    if (!truth_table || n < 0 || n > 24) return 0.0;
    int32_t N = 1 << n;
    double sum = 0.0;
    for (int32_t x = 0; x < N; x++) {
        int32_t parity = popcount64((uint64_t)x & subset) & 1;
        sum += (double)(truth_table[x] & 1) * (parity ? -1.0 : 1.0);
    }
    return sum / (double)N;
}

/* ============================================================================
 * Subset Convolution (Polynomial Multiplication / Function AND)
 *
 * The pointwise AND (f AND g)(x) = f(x) * g(x) corresponds to subset
 * convolution of their polynomial representations:
 *
 *   (f * g)_hat(S) = sum_{A union B = S, A cap B = empty} f_hat(A) * g_hat(B)
 *
 * This is the "subset convolution" operation. Using ranked Mobius transform,
 * it can be computed in O(n^2 * 2^n) time.
 *
 * Reference: Bjorklund, Husfeldt, Kaski, Koivisto (2007),
 *   "Fourier meets Mobius: fast subset convolution".
 *
 * Time: O(n^2 * 2^n). Practically limited to n <= 20.
 * ============================================================================ */

polynomial_t *poly_and(const polynomial_t *f, const polynomial_t *g) {
    if (!f || !g) return NULL;
    int32_t n = (f->num_vars > g->num_vars) ? f->num_vars : g->num_vars;
    int32_t N = 1 << n;
    if (N <= 0 || n > 20) return NULL;

    double *cf = (double *)calloc((size_t)N, sizeof(double));
    double *cg = (double *)calloc((size_t)N, sizeof(double));
    if (!cf || !cg) { free(cf); free(cg); return NULL; }

    for (int32_t i = 0; i < f->num_terms; i++) {
        uint64_t v = f->terms[i].mono.vars;
        if (v < (uint64_t)N) cf[v] += f->terms[i].coeff;
    }
    for (int32_t i = 0; i < g->num_terms; i++) {
        uint64_t v = g->terms[i].mono.vars;
        if (v < (uint64_t)N) cg[v] += g->terms[i].coeff;
    }

    /* Ranked transform approach */
    double **fr = (double **)calloc((size_t)(n + 1), sizeof(double *));
    double **gr = (double **)calloc((size_t)(n + 1), sizeof(double *));
    double **hr = (double **)calloc((size_t)(n + 1), sizeof(double *));
    for (int32_t k = 0; k <= n; k++) {
        fr[k] = (double *)calloc((size_t)N, sizeof(double));
        gr[k] = (double *)calloc((size_t)N, sizeof(double));
        hr[k] = (double *)calloc((size_t)N, sizeof(double));
        if (!fr[k] || !gr[k] || !hr[k]) {
            /* Cleanup on failure */
            for (int32_t kk = 0; kk <= n; kk++) {
                free(fr[kk]); free(gr[kk]); free(hr[kk]);
            }
            free(fr); free(gr); free(hr); free(cf); free(cg);
            return NULL;
        }
    }

    /* Split by rank */
    for (int32_t S = 0; S < N; S++) {
        int32_t rank = popcount64((uint64_t)S);
        fr[rank][S] = cf[S];
        gr[rank][S] = cg[S];
    }

    /* Zeta transform each rank */
    for (int32_t k = 0; k <= n; k++) {
        zeta_transform(fr[k], n);
        zeta_transform(gr[k], n);
    }

    /* Pointwise multiply ranks */
    for (int32_t S = 0; S < N; S++) {
        for (int32_t i = 0; i <= n; i++) {
            for (int32_t j = 0; j <= n - i; j++) {
                hr[i + j][S] += fr[i][S] * gr[j][S];
            }
        }
    }

    /* Mobius transform and extract correct-rank terms */
    for (int32_t k = 0; k <= n; k++) {
        mobius_transform(hr[k], n);
        for (int32_t S = 0; S < N; S++) {
            if (popcount64((uint64_t)S) == k) {
                cf[S] = hr[k][S];  /* reuse cf for result */
            }
        }
    }

    polynomial_t *result = poly_create(n, N, f->over_gf2 && g->over_gf2);
    if (!result) {
        for (int32_t k = 0; k <= n; k++) {
            free(fr[k]); free(gr[k]); free(hr[k]);
        }
        free(fr); free(gr); free(hr); free(cf); free(cg);
        return NULL;
    }
    for (int32_t S = 0; S < N; S++) {
        if (fabs(cf[S]) > 1e-12)
            poly_add_term(result, cf[S], (uint64_t)S);
    }

    for (int32_t k = 0; k <= n; k++) {
        free(fr[k]); free(gr[k]); free(hr[k]);
    }
    free(fr); free(gr); free(hr); free(cf); free(cg);
    return result;
}

/* ============================================================================
 * Polynomial XOR (GF(2) addition of Boolean functions)
 *
 * For GF(2) polynomials, addition IS XOR. So the polynomial
 * representing the XOR of two Boolean functions is simply
 * the sum of their individual polynomials over GF(2).
 * ============================================================================ */

polynomial_t *poly_xor(const polynomial_t *f, const polynomial_t *g) {
    return poly_add(f, g);
}

/* ============================================================================
 * Compute Polynomial Degree Vector (degree spectrum)
 *
 * For a Boolean function f, the degree spectrum shows how much "mass"
 * the function has at each degree level. This is useful for:
 * - Approximating f by low-degree polynomials
 * - Proving that f requires high-degree approximators
 * - Characterizing circuit complexity classes
 *
 * Returns array of size (n+1): d[0] = constant term mass,
 * d[1] = linear term mass, ..., d[n] = degree-n term mass.
 * ============================================================================ */

double *poly_degree_spectrum(const polynomial_t *p, int32_t n) {
    if (!p || n < 0) return NULL;
    double *spectrum = (double *)calloc((size_t)(n + 1), sizeof(double));
    if (!spectrum) return NULL;
    for (int32_t i = 0; i < p->num_terms; i++) {
        int32_t d = p->terms[i].mono.degree;
        if (d <= n) spectrum[d] += fabs(p->terms[i].coeff);
    }
    return spectrum;
}

/* ============================================================================
 * Polynomial Sparsity: count non-zero terms
 * ============================================================================ */

int32_t poly_sparsity(const polynomial_t *p) {
    if (!p) return 0;
    return p->num_terms;
}

/* ============================================================================
 * Check if polynomial is multi-linear (no variable squared)
 * In our bitmask representation this is always true, but we verify.
 * ============================================================================ */

bool poly_is_multilinear(const polynomial_t *p) {
    if (!p) return false;
    for (int32_t i = 0; i < p->num_terms; i++) {
        /* Our bitmask ensures each variable appears at most once */
        /* But check: degree = popcount(var_mask) */
        if (popcount64(p->terms[i].mono.vars) != p->terms[i].mono.degree)
            return false;
    }
    return true;
}
