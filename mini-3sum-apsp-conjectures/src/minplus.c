/* minplus.c - Min-Plus (Tropical) Matrix Multiplication */
#include "minplus.h"
#include "apsp.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

static double inf_min(double a, double b) { return a < b ? a : b; }

mp_matrix_t *mp_matrix_create(int32_t rows, int32_t cols) {
    mp_matrix_t *m = malloc(sizeof(mp_matrix_t));
    if (!m) return NULL;
    m->rows = rows; m->cols = cols;
    m->data = malloc((size_t)(rows * cols) * sizeof(double));
    if (!m->data) { free(m); return NULL; }
    for (int32_t i = 0; i < rows * cols; i++) m->data[i] = MP_INF;
    return m;
}

void mp_matrix_destroy(mp_matrix_t *m) {
    if (!m) return;
    free(m->data); free(m);
}

double mp_get(const mp_matrix_t *m, int32_t i, int32_t j) {
    if (!m || i < 0 || i >= m->rows || j < 0 || j >= m->cols) return MP_INF;
    return m->data[i * m->cols + j];
}

void mp_set(mp_matrix_t *m, int32_t i, int32_t j, double val) {
    if (!m || i < 0 || i >= m->rows || j < 0 || j >= m->cols) return;
    m->data[i * m->cols + j] = val;
}

void mp_fill(mp_matrix_t *m, double val) {
    if (!m) return;
    for (int32_t i = 0; i < m->rows * m->cols; i++) m->data[i] = val;
}

mp_matrix_t *mp_naive_multiply(const mp_matrix_t *A, const mp_matrix_t *B) {
    if (!A || !B || A->cols != B->rows) return NULL;
    int32_t n = A->rows, m = B->cols, p = A->cols;
    mp_matrix_t *C = mp_matrix_create(n, m);
    if (!C) return NULL;
    for (int32_t i = 0; i < n; i++) {
        for (int32_t j = 0; j < m; j++) {
            double best = MP_INF;
            for (int32_t k = 0; k < p; k++) {
                double aik = A->data[i * p + k];
                double bkj = B->data[k * m + j];
                if (aik < MP_INF / 2.0 && bkj < MP_INF / 2.0)
                    best = inf_min(best, aik + bkj);
            }
            C->data[i * m + j] = best;
        }
    }
    return C;
}

mp_matrix_t *mp_strassen_inspired(const mp_matrix_t *A, const mp_matrix_t *B) {
    if (!A || !B || A->rows != A->cols || B->rows != B->cols || A->rows != B->rows)
        return mp_naive_multiply(A, B);
    int32_t n = A->rows;
    if (n <= 64) return mp_naive_multiply(A, B);
    int32_t half = n / 2;
    mp_matrix_t *C = mp_matrix_create(n, n);
    if (!C) return NULL;
    mp_matrix_t *A11 = mp_matrix_create(half, half);
    mp_matrix_t *A12 = mp_matrix_create(half, half);
    mp_matrix_t *A21 = mp_matrix_create(half, half);
    mp_matrix_t *A22 = mp_matrix_create(half, half);
    mp_matrix_t *B11 = mp_matrix_create(half, half);
    mp_matrix_t *B12 = mp_matrix_create(half, half);
    mp_matrix_t *B21 = mp_matrix_create(half, half);
    mp_matrix_t *B22 = mp_matrix_create(half, half);
    for (int32_t i = 0; i < half; i++) {
        for (int32_t j = 0; j < half; j++) {
            A11->data[i*half+j] = A->data[i*n+j];
            A12->data[i*half+j] = A->data[i*n+j+half];
            A21->data[i*half+j] = A->data[(i+half)*n+j];
            A22->data[i*half+j] = A->data[(i+half)*n+j+half];
            B11->data[i*half+j] = B->data[i*n+j];
            B12->data[i*half+j] = B->data[i*n+j+half];
            B21->data[i*half+j] = B->data[(i+half)*n+j];
            B22->data[i*half+j] = B->data[(i+half)*n+j+half];
        }
    }
    mp_matrix_t *M1 = mp_strassen_inspired(A11, B11);
    mp_matrix_t *M2 = mp_strassen_inspired(A12, B21);
    mp_matrix_t *M3 = mp_strassen_inspired(A11, B12);
    mp_matrix_t *M4 = mp_strassen_inspired(A12, B22);
    mp_matrix_t *M5 = mp_strassen_inspired(A21, B11);
    mp_matrix_t *M6 = mp_strassen_inspired(A22, B21);
    mp_matrix_t *M7 = mp_strassen_inspired(A21, B12);
    mp_matrix_t *M8 = mp_strassen_inspired(A22, B22);
    for (int32_t i = 0; i < half; i++) {
        for (int32_t j = 0; j < half; j++) {
            C->data[i*n+j] = inf_min(M1->data[i*half+j], M2->data[i*half+j]);
            C->data[i*n+j+half] = inf_min(M3->data[i*half+j], M4->data[i*half+j]);
            C->data[(i+half)*n+j] = inf_min(M5->data[i*half+j], M6->data[i*half+j]);
            C->data[(i+half)*n+j+half] = inf_min(M7->data[i*half+j], M8->data[i*half+j]);
        }
    }
    mp_matrix_destroy(A11); mp_matrix_destroy(A12);
    mp_matrix_destroy(A21); mp_matrix_destroy(A22);
    mp_matrix_destroy(B11); mp_matrix_destroy(B12);
    mp_matrix_destroy(B21); mp_matrix_destroy(B22);
    mp_matrix_destroy(M1); mp_matrix_destroy(M2); mp_matrix_destroy(M3); mp_matrix_destroy(M4);
    mp_matrix_destroy(M5); mp_matrix_destroy(M6); mp_matrix_destroy(M7); mp_matrix_destroy(M8);
    return C;
}

mp_matrix_t *mp_repeated_squaring(const mp_matrix_t *A, int32_t k) {
    if (!A || A->rows != A->cols || k < 1) return NULL;
    int32_t n = A->rows;
    mp_matrix_t *result = mp_matrix_create(n, n);
    if (!result) return NULL;
    for (int32_t i = 0; i < n; i++) result->data[i*n+i] = 0.0;
    mp_matrix_t *base = mp_naive_multiply(A, A);
    if (!base) { mp_matrix_destroy(result); return NULL; }
    mp_matrix_t *tmp = base;
    for (int32_t i = 0; i < k; i++) {
        mp_matrix_t *next = mp_naive_multiply(result, tmp);
        mp_matrix_destroy(result);
        result = next;
        if (!result) break;
    }
    mp_matrix_destroy(base);
    return result;
}

mp_matrix_t *mp_fast_apsp(const mp_matrix_t *A) {
    if (!A || A->rows != A->cols) return NULL;
    int32_t n = A->rows;
    mp_matrix_t *D = mp_matrix_create(n, n);
    if (!D) return NULL;
    for (int32_t i = 0; i < n * n; i++) D->data[i] = A->data[i];
    for (int32_t k = 0; k < n; k++) {
        for (int32_t i = 0; i < n; i++) {
            if (D->data[i*n+k] >= MP_INF/2.0) continue;
            for (int32_t j = 0; j < n; j++) {
                if (D->data[k*n+j] >= MP_INF/2.0) continue;
                double nd = D->data[i*n+k] + D->data[k*n+j];
                if (nd < D->data[i*n+j]) D->data[i*n+j] = nd;
            }
        }
    }
    return D;
}

mp_matrix_t *mp_closure(const mp_matrix_t *A) {
    return mp_fast_apsp(A);
}

bool mp_is_idempotent(const mp_matrix_t *A) {
    if (!A || A->rows != A->cols) return false;
    mp_matrix_t *A2 = mp_naive_multiply(A, A);
    if (!A2) return false;
    bool result = true;
    for (int32_t i = 0; i < A->rows && result; i++)
        for (int32_t j = 0; j < A->cols && result; j++)
            if (fabs(A->data[i*A->cols+j] - A2->data[i*A2->cols+j]) > 1e-12)
                result = false;
    mp_matrix_destroy(A2);
    return result;
}

void mp_print(const mp_matrix_t *m) {
    if (!m) { printf("NULL\n"); return; }
    printf("Min-plus matrix %dx%d:\n", m->rows, m->cols);
    int32_t dr = m->rows < 10 ? m->rows : 10;
    int32_t dc = m->cols < 10 ? m->cols : 10;
    for (int32_t i = 0; i < dr; i++) {
        for (int32_t j = 0; j < dc; j++) {
            double v = m->data[i * m->cols + j];
            if (v >= MP_INF / 2.0) printf("  INF");
            else printf(" %4.1f", v);
        }
        printf("\n");
    }
    if (m->rows > 10 || m->cols > 10) printf("  ... (truncated)\n");
}

/* === Additional Min-Plus Operations ================================== */

mp_matrix_t *mp_matrix_from_adjacency(const apsp_adjacency_t *adj) {
    if (!adj) return NULL;
    mp_matrix_t *m = mp_matrix_create(adj->n, adj->n);
    if (!m) return NULL;
    for (int32_t i = 0; i < adj->n * adj->n; i++)
        m->data[i] = (adj->weights[i] >= APSP_INF / 2.0) ? MP_INF : adj->weights[i];
    return m;
}

mp_matrix_t *mp_add(const mp_matrix_t *A, const mp_matrix_t *B) {
    if (!A || !B || A->rows != B->rows || A->cols != B->cols) return NULL;
    mp_matrix_t *C = mp_matrix_create(A->rows, A->cols);
    if (!C) return NULL;
    for (int32_t i = 0; i < A->rows * A->cols; i++)
        C->data[i] = (A->data[i] < B->data[i]) ? A->data[i] : B->data[i];
    return C;
}

mp_matrix_t *mp_transpose(const mp_matrix_t *A) {
    if (!A) return NULL;
    mp_matrix_t *T = mp_matrix_create(A->cols, A->rows);
    if (!T) return NULL;
    for (int32_t i = 0; i < A->rows; i++)
        for (int32_t j = 0; j < A->cols; j++)
            T->data[j * T->cols + i] = A->data[i * A->cols + j];
    return T;
}

bool mp_is_symmetric(const mp_matrix_t *A) {
    if (!A || A->rows != A->cols) return false;
    for (int32_t i = 0; i < A->rows; i++)
        for (int32_t j = i + 1; j < A->cols; j++)
            if (fabs(A->data[i * A->cols + j] - A->data[j * A->cols + i]) > 1e-12)
                return false;
    return true;
}

mp_matrix_t *mp_identity(int32_t n) {
    mp_matrix_t *I = mp_matrix_create(n, n);
    if (!I) return NULL;
    for (int32_t i = 0; i < n; i++) I->data[i * n + i] = 0.0;
    return I;
}

bool mp_equal(const mp_matrix_t *A, const mp_matrix_t *B) {
    if (!A || !B) return A == B;
    if (A->rows != B->rows || A->cols != B->cols) return false;
    for (int32_t i = 0; i < A->rows * A->cols; i++)
        if (fabs(A->data[i] - B->data[i]) > 1e-12) return false;
    return true;
}

double mp_min_entry(const mp_matrix_t *A) {
    if (!A || A->rows == 0 || A->cols == 0) return MP_INF;
    double mv = A->data[0];
    for (int32_t i = 1; i < A->rows * A->cols; i++)
        if (A->data[i] < mv) mv = A->data[i];
    return mv;
}

double mp_max_finite_entry(const mp_matrix_t *A) {
    if (!A || A->rows == 0 || A->cols == 0) return 0.0;
    double mv = -MP_INF;
    bool found = false;
    for (int32_t i = 0; i < A->rows * A->cols; i++) {
        if (A->data[i] < MP_INF / 2.0) {
            if (A->data[i] > mv) mv = A->data[i];
            found = true;
        }
    }
    return found ? mv : 0.0;
}

mp_matrix_t *mp_closure_via_floyd(const mp_matrix_t *A) {
    return mp_fast_apsp(A);
}

int32_t mp_count_finite(const mp_matrix_t *A) {
    if (!A) return 0;
    int32_t count = 0;
    for (int32_t i = 0; i < A->rows * A->cols; i++)
        if (A->data[i] < MP_INF / 2.0) count++;
    return count;
}

double mp_density(const mp_matrix_t *A) {
    if (!A || A->rows == 0 || A->cols == 0) return 0.0;
    return (double)mp_count_finite(A) / (A->rows * A->cols);
}

mp_matrix_t *mp_copy(const mp_matrix_t *A) {
    if (!A) return NULL;
    mp_matrix_t *C = mp_matrix_create(A->rows, A->cols);
    if (!C) return NULL;
    memcpy(C->data, A->data, (size_t)(A->rows * A->cols) * sizeof(double));
    return C;
}

mp_matrix_t *mp_shortest_path_via_minplus(mp_matrix_t *A, int32_t src, int32_t dst) {
    if (!A || A->rows != A->cols) return NULL;
    int32_t n = A->rows;
    mp_matrix_t *D = mp_fast_apsp(A);
    if (!D) return NULL;
    double dist = D->data[src * n + dst];
    mp_matrix_t *result = mp_matrix_create(1, 1);
    if (result) result->data[0] = dist;
    mp_matrix_destroy(D);
    return result;
}

mp_matrix_t *mp_from_array(const double *data, int32_t n) {
    mp_matrix_t *m = mp_matrix_create(n, n);
    if (!m) return NULL;
    memcpy(m->data, data, (size_t)(n * n) * sizeof(double));
    return m;
}

void mp_export_to_array(const mp_matrix_t *m, double *out) {
    if (!m || !out) return;
    memcpy(out, m->data, (size_t)(m->rows * m->cols) * sizeof(double));
}

mp_matrix_t *mp_elementwise_min(const mp_matrix_t *A, const mp_matrix_t *B) {
    if (!A || !B || A->rows != B->rows || A->cols != B->cols) return NULL;
    mp_matrix_t *C = mp_matrix_create(A->rows, A->cols);
    if (!C) return NULL;
    for (int32_t i = 0; i < A->rows * A->cols; i++)
        C->data[i] = (A->data[i] < B->data[i]) ? A->data[i] : B->data[i];
    return C;
}

mp_matrix_t *mp_elementwise_add(const mp_matrix_t *A, const mp_matrix_t *B) {
    if (!A || !B || A->rows != B->rows || A->cols != B->cols) return NULL;
    mp_matrix_t *C = mp_matrix_create(A->rows, A->cols);
    if (!C) return NULL;
    for (int32_t i = 0; i < A->rows * A->cols; i++) {
        double a = A->data[i], b = B->data[i];
        if (a >= MP_INF / 2.0 || b >= MP_INF / 2.0)
            C->data[i] = MP_INF;
        else
            C->data[i] = a + b;
    }
    return C;
}

void mp_scale(mp_matrix_t *A, double factor) {
    if (!A) return;
    for (int32_t i = 0; i < A->rows * A->cols; i++)
        if (A->data[i] < MP_INF / 2.0) A->data[i] *= factor;
}

mp_matrix_t *mp_random(int32_t n, double density, double max_weight, uint64_t seed) {
    if (seed == 0) seed = (uint64_t)time(NULL);
    srand((unsigned int)seed);
    mp_matrix_t *m = mp_matrix_create(n, n);
    if (!m) return NULL;
    for (int32_t i = 0; i < n; i++) {
        m->data[i * n + i] = 0.0;
        for (int32_t j = 0; j < n; j++) {
            if (i != j && ((double)rand() / RAND_MAX) < density)
                m->data[i * n + j] = ((double)rand() / RAND_MAX) * max_weight + 0.1;
        }
    }
    return m;
}

mp_matrix_t *mp_submatrix(const mp_matrix_t *A, int32_t r0, int32_t r1, int32_t c0, int32_t c1) {
    if (!A || r0<0 || r1> A->rows || c0<0 || c1> A->cols || r0>=r1 || c0>=c1) return NULL;
    int32_t nr=r1-r0, nc=c1-c0;
    mp_matrix_t *S=mp_matrix_create(nr,nc);
    if (!S) return NULL;
    for (int32_t i=0; i<nr; i++)
        for (int32_t j=0; j<nc; j++)
            S->data[i*nc+j]=A->data[(r0+i)*A->cols+(c0+j)];
    return S;
}
void mp_set_submatrix(mp_matrix_t *A, int32_t r0, int32_t c0, const mp_matrix_t *B) {
    if (!A || !B || r0<0 || c0<0 || r0+B->rows>A->rows || c0+B->cols>A->cols) return;
    for (int32_t i=0; i<B->rows; i++)
        for (int32_t j=0; j<B->cols; j++)
            A->data[(r0+i)*A->cols+(c0+j)]=B->data[i*B->cols+j];
}
bool mp_is_triangle_inequality(const mp_matrix_t *A) {
    if (!A || A->rows!=A->cols) return false;
    for (int32_t i=0; i<A->rows; i++)
        for (int32_t j=0; j<A->cols; j++)
            for (int32_t k=0; k<A->cols; k++)
                if (A->data[i*A->cols+k]<MP_INF/2.0 && A->data[k*A->cols+j]<MP_INF/2.0)
                    if (A->data[i*A->cols+j] > A->data[i*A->cols+k]+A->data[k*A->cols+j]+1e-12)
                        return false;
    return true;
}

bool mp_verify_apsp_equivalence(const mp_matrix_t *A) {
    if (!A || A->rows != A->cols) return false;
    mp_matrix_t *D = mp_fast_apsp(A);
    if (!D) return false;
    mp_matrix_t *D2 = mp_naive_multiply(D, D);
    if (!D2) { mp_matrix_destroy(D); return false; }
    bool ok = mp_equal(D, D2);
    printf("APSP equivalence check (D = D*D): %s\n", ok ? "PASS" : "FAIL");
    mp_matrix_destroy(D); mp_matrix_destroy(D2);
    return ok;
}
void mp_self_test(void) {
    printf("=== Min-Plus Self-Test ===\n");
    int32_t n = 4;
    mp_matrix_t *A = mp_random(n, 0.5, 10.0, 42);
    mp_matrix_t *B = mp_naive_multiply(A, A);
    printf("A (4x4 random, density=0.5):\n"); mp_print(A);
    printf("A*A (min-plus product):\n"); mp_print(B);
    mp_matrix_t *D = mp_fast_apsp(A);
    printf("APSP via Floyd (closure):\n"); mp_print(D);
    printf("Is idempotent: %s\n", mp_is_idempotent(D) ? "yes" : "no");
    mp_matrix_destroy(A); mp_matrix_destroy(B); mp_matrix_destroy(D);
}
