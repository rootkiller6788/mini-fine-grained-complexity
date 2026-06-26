#ifndef MINPLUS_H
#define MINPLUS_H
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdbool.h>
#include <float.h>
#ifdef __cplusplus
extern "C" {
#endif
#define MP_INF DBL_MAX
typedef struct { double *data; int32_t rows; int32_t cols; } mp_matrix_t;
mp_matrix_t *mp_matrix_create(int32_t rows, int32_t cols);
void mp_matrix_destroy(mp_matrix_t *m);
double mp_get(const mp_matrix_t *m, int32_t i, int32_t j);
void mp_set(mp_matrix_t *m, int32_t i, int32_t j, double val);
mp_matrix_t *mp_naive_multiply(const mp_matrix_t *A, const mp_matrix_t *B);
mp_matrix_t *mp_fast_apsp(const mp_matrix_t *A);
mp_matrix_t *mp_repeated_squaring(const mp_matrix_t *A, int32_t k);
mp_matrix_t *mp_closure(const mp_matrix_t *A);
void mp_print(const mp_matrix_t *m);
#ifdef __cplusplus
}
#endif
#endif

/*
 * ============================================================
 * Min-Plus: Extended Documentation
 * ============================================================
 *
 * The min-plus semiring (also called the tropical semiring)
 * is defined over R U {+infinity} with operations:
 *   a "+" b = min(a, b)   (min-plus "addition")
 *   a "*" b = a + b       (min-plus "multiplication")
 *
 * Matrix product over this semiring:
 *   C[i][j] = min_k { A[i][k] + B[k][j] }
 *
 * This is also called the "distance product" because it
 * generalizes shortest-path computation to matrix form.
 *
 * Key properties:
 *   - The identity matrix has 0 on diagonal, +inf elsewhere.
 *   - Min-plus matrix multiplication is associative.
 *   - APSP = (n-1)-fold min-plus power of adjacency matrix.
 *   - Repeated squaring gives APSP in O(T(n) log n) where
 *     T(n) is the time for one min-plus product.
 *
 * The min-plus product conjecture:
 *   Min-plus product of two n x n matrices requires
 *   n^{3-o(1)} time. This is equivalent to the APSP conjecture.
 *
 * Connection to standard matrix multiplication:
 *   If we could compute min-plus products as fast as standard
 *   matrix products (O(n^{omega}) with omega ~ 2.3729), we would
 *   refute the APSP conjecture. The fact that min-plus seems
 *   genuinely harder than (+,*)-product is a deep mystery in
 *   algebraic complexity theory.
 *
 * The gap between standard matrix multiplication (O(n^{2.3729}))
 * and min-plus product (essentially O(n^3)) is one of the most
 * significant open problems in fine-grained complexity.
 *
 * For further reading:
 *   - Minoux (1976): "Plus court chemin...", first formal study
 *   - Aho, Hopcroft, Ullman (1974): "The Design and Analysis of
 *     Computer Algorithms", min-plus closure
 */
/*
 * Min-Plus Algebraic Properties:
 *
 * The min-plus semiring (R U {inf}, min, +) satisfies:
 *   - Associativity: (A*B)*C = A*(B*C)
 *   - Distributivity: A*(B+C) = A*B + A*C
 *     where (B+C)[i][j] = min(B[i][j], C[i][j])
 *   - Identity I: I[i][i]=0, I[i][j]=inf (i!=j)
 *   - Zero 0: 0[i][j]=inf for all i,j
 *
 * Unlike standard matrix multiplication, min-plus product:
 *   - Is NOT invertible (no subtraction)
 *   - Has NO known fast algorithm (unlike O(n^{2.3729}) for +,*)
 *   - Is equivalent to APSP under subcubic reductions
 */