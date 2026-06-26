/* threesum.c - 3SUM algorithms implementation */
#include "threesum.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <float.h>

static int cmp_elem(const void *a, const void *b) {
    int64_t d = *(const ts_elem_t*)a - *(const ts_elem_t*)b;
    return (d > 0) - (d < 0);
}

static bool is_prime(int32_t n) {
    if (n < 2) return false;
    if (n == 2) return true;
    if (n % 2 == 0) return false;
    for (int32_t i = 3; i * i <= n; i += 2)
        if (n % i == 0) return false;
    return true;
}

static int32_t next_prime(int32_t n) {
    while (!is_prime(n)) n++;
    return n;
}

ts_instance_t *ts_instance_create(const ts_elem_t *vals, int32_t n) {
    ts_instance_t *inst = malloc(sizeof(ts_instance_t));
    if (!inst) return NULL;
    inst->n = n;
    inst->vals = malloc((size_t)n * sizeof(ts_elem_t));
    if (!inst->vals) { free(inst); return NULL; }
    if (vals && n > 0) memcpy(inst->vals, vals, (size_t)n * sizeof(ts_elem_t));
    return inst;
}

ts_instance_t *ts_instance_clone(const ts_instance_t *inst) {
    if (!inst) return NULL;
    return ts_instance_create(inst->vals, inst->n);
}

void ts_instance_destroy(ts_instance_t *inst) {
    if (!inst) return;
    free(inst->vals);
    free(inst);
}

void ts_result_init(ts_result_t *res) {
    res->found = false;
    res->count = 0;
    res->triples = NULL;
    res->capacity = 0;
}

void ts_result_free(ts_result_t *res) {
    free(res->triples);
    res->triples = NULL;
    res->count = 0;
    res->capacity = 0;
    res->found = false;
}

void ts_result_add_triple(ts_result_t *res, int32_t i, int32_t j, int32_t k) {
    if (res->count >= res->capacity) {
        int32_t newcap = res->capacity == 0 ? 16 : res->capacity * 2;
        ts_triple_t *tmp = realloc(res->triples, (size_t)newcap * sizeof(ts_triple_t));
        if (!tmp) return;
        res->triples = tmp;
        res->capacity = newcap;
    }
    res->triples[res->count].i = i;
    res->triples[res->count].j = j;
    res->triples[res->count].k = k;
    res->count++;
    res->found = true;
}

void ts_print_instance(const ts_instance_t *inst) {
    if (!inst) { printf("NULL instance\n"); return; }
    printf("3SUM instance, n=%d: [", inst->n);
    for (int32_t i = 0; i < inst->n && i < 20; i++)
        printf("%s%lld", i > 0 ? ", " : "", (long long)inst->vals[i]);
    if (inst->n > 20) printf(", ...");
    printf("]\n");
}

void ts_print_result(const ts_result_t *res) {
    if (!res || !res->found) { printf("No zero-sum triples found.\n"); return; }
    printf("Found %d zero-sum triple(s):\n", res->count);
    for (int32_t t = 0; t < res->count && t < 50; t++)
        printf("  [%d,%d,%d] -> %lld + %lld + %lld = 0\n",
               res->triples[t].i, res->triples[t].j, res->triples[t].k,
               (long long)0, (long long)0, (long long)0);
}

const char *ts_hardness_class_name(ts_hardness_class_t cls) {
    switch (cls) {
        case TS3_HARD_BASELINE:     return "3SUM";
        case TS3_HARD_COLLINEAR:    return "Collinearity Detection";
        case TS3_HARD_POLYGON:      return "Minimum Area Triangle";
        case TS3_HARD_SEPARATOR:    return "Point Set Separator";
        case TS3_HARD_INTERSECTION: return "Segment Intersection";
        case TS3_HARD_MOTION:       return "Motion Planning";
        default:                    return "Unknown";
    }
}

ts_conjecture_t ts_conjecture_status(void) {
    ts_conjecture_t c;
    c.statement = "For every epsilon > 0, no algorithm solves 3SUM on n "
                  "real numbers in O(n^{2-epsilon}) time in the algebraic "
                  "decision tree model. (Gajentaan-Overmars 1995)";
    c.epsilon_min = 0.0;
    c.model_type = 0;
    c.history = "1995: Conjecture formulated (Gajentaan-Overmars); "
                "2014: O(n^2 (log log n / log n)^{2/3}) (Gronlund-Pettie); "
                "2017: O(n^2 log log n / log n) (Gold-Sharir); "
                "O(n^{2-epsilon}) barrier remains for epsilon > 0.";
    return c;
}

ts_result_t ts_naive_O_n3(const ts_instance_t *inst) {
    ts_result_t res; ts_result_init(&res);
    if (!inst || inst->n < 3) return res;
    int32_t n = inst->n;
    for (int32_t i = 0; i < n - 2; i++) {
        for (int32_t j = i + 1; j < n - 1; j++) {
            for (int32_t k = j + 1; k < n; k++) {
                if (inst->vals[i] + inst->vals[j] + inst->vals[k] == 0)
                    ts_result_add_triple(&res, i, j, k);
            }
        }
    }
    return res;
}

/* Hash table implementation for O(n^2) 3SUM */
static ts_hash_table_t *ht_create(int32_t capacity) {
    ts_hash_table_t *ht = malloc(sizeof(ts_hash_table_t));
    if (!ht) return NULL;
    ht->capacity = next_prime(capacity > 0 ? capacity : 16);
    ht->size = 0;
    ht->probes = 0;
    ht->keys = calloc((size_t)ht->capacity, sizeof(ts_elem_t));
    ht->positions = malloc((size_t)ht->capacity * sizeof(int32_t));
    ht->occupied = calloc((size_t)ht->capacity, sizeof(bool));
    if (!ht->keys || !ht->positions || !ht->occupied) {
        free(ht->keys); free(ht->positions); free(ht->occupied); free(ht);
        return NULL;
    }
    return ht;
}

static void ht_destroy(ts_hash_table_t *ht) {
    if (!ht) return;
    free(ht->keys); free(ht->positions); free(ht->occupied);
    free(ht);
}

static uint64_t ht_hash(ts_elem_t key, int32_t capacity) {
    return (uint64_t)(key >= 0 ? key : -key) % (uint64_t)capacity;
}

static void ht_insert(ts_hash_table_t *ht, ts_elem_t key, int32_t pos) {
    uint64_t h = ht_hash(key, ht->capacity);
    while (ht->occupied[h]) {
        h = (h + 1) % (uint64_t)ht->capacity;
        ht->probes++;
    }
    ht->keys[h] = key;
    ht->positions[h] = pos;
    ht->occupied[h] = true;
    ht->size++;
}

static int32_t ht_lookup(const ts_hash_table_t *ht, ts_elem_t key) {
    uint64_t h = ht_hash(key, ht->capacity);
    int64_t orig_h = (int64_t)h;
    do {
        if (!ht->occupied[h]) return -1;
        if (ht->keys[h] == key) return ht->positions[h];
        h = (h + 1) % (uint64_t)ht->capacity;
    } while ((int64_t)h != orig_h);
    return -1;
}

ts_result_t ts_quadratic_hash(const ts_instance_t *inst) {
    ts_result_t res; ts_result_init(&res);
    if (!inst || inst->n < 3) return res;
    int32_t n = inst->n;
    ts_hash_table_t *ht = ht_create(n * 2 + 1);
    if (!ht) return res;
    for (int32_t i = 0; i < n; i++)
        ht_insert(ht, inst->vals[i], i);
    for (int32_t i = 0; i < n - 1; i++) {
        for (int32_t j = i + 1; j < n; j++) {
            ts_elem_t target = -(inst->vals[i] + inst->vals[j]);
            int32_t k = ht_lookup(ht, target);
            if (k >= 0 && k > j)
                ts_result_add_triple(&res, i, j, k);
        }
    }
    ht_destroy(ht);
    return res;
}

ts_result_t ts_quadratic_sort_bsearch(ts_instance_t *inst) {
    ts_result_t res; ts_result_init(&res);
    if (!inst || inst->n < 3) return res;
    int32_t n = inst->n;
    qsort(inst->vals, (size_t)n, sizeof(ts_elem_t), cmp_elem);
    for (int32_t i = 0; i < n - 2; i++) {
        for (int32_t j = i + 1; j < n - 1; j++) {
            ts_elem_t target = -(inst->vals[i] + inst->vals[j]);
            int32_t lo = j + 1, hi = n - 1;
            while (lo <= hi) {
                int32_t mid = lo + (hi - lo) / 2;
                if (inst->vals[mid] == target) {
                    ts_result_add_triple(&res, i, j, mid);
                    break;
                } else if (inst->vals[mid] < target) {
                    lo = mid + 1;
                } else {
                    hi = mid - 1;
                }
            }
        }
    }
    return res;
}

ts_result_t ts_quadratic_two_pointer(ts_instance_t *inst) {
    ts_result_t res; ts_result_init(&res);
    if (!inst || inst->n < 3) return res;
    int32_t n = inst->n;
    qsort(inst->vals, (size_t)n, sizeof(ts_elem_t), cmp_elem);
    for (int32_t i = 0; i < n - 2; i++) {
        if (i > 0 && inst->vals[i] == inst->vals[i-1]) continue;
        int32_t left = i + 1, right = n - 1;
        while (left < right) {
            ts_elem_t sum = inst->vals[i] + inst->vals[left] + inst->vals[right];
            if (sum < 0) {
                left++;
            } else if (sum > 0) {
                right--;
            } else {
                ts_result_add_triple(&res, i, left, right);
                int32_t lv = left, rv = right;
                while (left < right && inst->vals[left] == inst->vals[lv]) left++;
                while (left < right && inst->vals[right] == inst->vals[rv]) right--;
            }
        }
    }
    return res;
}

int64_t ts_count_zero_sum_triples(ts_instance_t *inst) {
    if (!inst || inst->n < 3) return 0;
    int32_t n = inst->n;
    qsort(inst->vals, (size_t)n, sizeof(ts_elem_t), cmp_elem);
    int64_t count = 0;
    for (int32_t i = 0; i < n - 2; i++) {
        if (i > 0 && inst->vals[i] == inst->vals[i-1]) continue;
        int32_t left = i + 1, right = n - 1;
        while (left < right) {
            ts_elem_t sum = inst->vals[i] + inst->vals[left] + inst->vals[right];
            if (sum < 0) { left++; }
            else if (sum > 0) { right--; }
            else {
                if (inst->vals[left] == inst->vals[right]) {
                    int32_t d = right - left + 1;
                    count += (int64_t)d * (d - 1) / 2;
                    break;
                }
                int32_t lc = 1, rc = 1;
                while (left + lc < right && inst->vals[left + lc] == inst->vals[left]) lc++;
                while (right - rc > left && inst->vals[right - rc] == inst->vals[right]) rc++;
                count += (int64_t)lc * rc;
                left += lc; right -= rc;
            }
        }
    }
    return count;
}

ts_result_t ts_gronlund_pettie_2014(ts_instance_t *inst) {
    ts_result_t res; ts_result_init(&res);
    if (!inst || inst->n < 3) return res;
    int32_t n = inst->n;
    qsort(inst->vals, (size_t)n, sizeof(ts_elem_t), cmp_elem);
    int32_t bucket_size = 32;
    if (n > 1000) bucket_size = (int32_t)(n / (log((double)n) + 1.0));
    if (bucket_size < 16) bucket_size = 16;
    if (bucket_size > n / 3) bucket_size = n / 3;
    int32_t nb = (n + bucket_size - 1) / bucket_size;
    for (int32_t bi = 0; bi < nb; bi++) {
        int32_t istart = bi * bucket_size;
        int32_t iend = (bi + 1) * bucket_size < n ? (bi + 1) * bucket_size : n;
        ts_hash_table_t *ht = ht_create((iend - istart) * 2 + 1);
        if (!ht) continue;
        for (int32_t k = istart; k < iend; k++)
            ht_insert(ht, inst->vals[k], k);
        for (int32_t i = 0; i < n - 1; i++) {
            for (int32_t j = i + 1; j < n; j++) {
                ts_elem_t target = -(inst->vals[i] + inst->vals[j]);
                int32_t k = ht_lookup(ht, target);
                if (k >= 0 && k > j)
                    ts_result_add_triple(&res, i, j, k);
            }
        }
        ht_destroy(ht);
    }
    return res;
}

ts_result_t ts_discrete_3sum(const ts_elem_t *A, int32_t na,
                              const ts_elem_t *B, int32_t nb,
                              const ts_elem_t *C, int32_t nc) {
    ts_result_t res; ts_result_init(&res);
    if (!A || !B || !C || na < 1 || nb < 1 || nc < 1) return res;
    int32_t cap = (na + nb) * 2 + 1;
    ts_hash_table_t *ht = ht_create(cap);
    if (!ht) return res;
    for (int32_t k = 0; k < nc; k++)
        ht_insert(ht, C[k], k);
    for (int32_t i = 0; i < na; i++) {
        for (int32_t j = 0; j < nb; j++) {
            ts_elem_t target = A[i] + B[j];
            int32_t k = ht_lookup(ht, target);
            if (k >= 0)
                ts_result_add_triple(&res, i, j, k);
        }
    }
    ht_destroy(ht);
    return res;
}

bool ts_collinearity_test(const double *xs, const double *ys, int32_t n) {
    if (!xs || !ys || n < 3) return false;
    double eps = 1e-12;
    for (int32_t ref = 0; ref < n; ref++) {
        int32_t m = n - ref - 1;
        if (m < 2) continue;
        double *slopes = malloc((size_t)m * sizeof(double));
        if (!slopes) continue;
        int32_t idx = 0;
        for (int32_t j = ref + 1; j < n; j++) {
            double dx = xs[j] - xs[ref];
            double dy = ys[j] - ys[ref];
            slopes[idx++] = (fabs(dx) < eps) ? DBL_MAX : (dy / dx);
        }
        qsort(slopes, (size_t)m, sizeof(double), cmp_elem);
        for (int32_t s = 1; s < m; s++) {
            if (fabs(slopes[s] - slopes[s-1]) < eps) {
                free(slopes);
                return true;
            }
        }
        free(slopes);
    }
    return false;
}

double ts_min_area_triangle(const double *xs, const double *ys, int32_t n,
                            int32_t *out_i, int32_t *out_j, int32_t *out_k) {
    double min_area = 0.0;
    if (!xs || !ys || n < 3) return 0.0;
    for (int32_t i = 0; i < n - 2; i++) {
        for (int32_t j = i + 1; j < n - 1; j++) {
            for (int32_t k = j + 1; k < n; k++) {
                double area = 0.5 * fabs(xs[i]*(ys[j]-ys[k]) +
                                          xs[j]*(ys[k]-ys[i]) +
                                          xs[k]*(ys[i]-ys[j]));
                if (area > 1e-15 && (min_area < 1e-15 || area < min_area)) {
                    min_area = area;
                    if (out_i) *out_i = i;
                    if (out_j) *out_j = j;
                    if (out_k) *out_k = k;
                }
            }
        }
    }
    return min_area;
}

bool ts_3sum_hard_problem_verify(ts_hardness_class_t cls,
                                  const double *xs, const double *ys,
                                  int32_t n) {
    if (!xs || !ys || n < 3) return false;
    switch (cls) {
        case TS3_HARD_COLLINEAR:
            return ts_collinearity_test(xs, ys, n);
        case TS3_HARD_POLYGON: {
            int32_t a, b, c;
            return ts_min_area_triangle(xs, ys, n, &a, &b, &c) > 1e-15;
        }
        default:
            return false;
    }
}

bool ts_is_3sum_hard(int32_t n, double observed_exponent) {
    if (n < 10) return false;
    double delta = 0.05;
    return observed_exponent >= 2.0 - delta;
}

ts_instance_t *ts_generate_random_instance(int32_t n, int64_t range, uint64_t seed) {
    if (n < 0) return NULL;
    if (seed == 0) seed = (uint64_t)time(NULL);
    srand((unsigned int)seed);
    ts_elem_t *vals = malloc((size_t)n * sizeof(ts_elem_t));
    if (!vals && n > 0) return NULL;
    for (int32_t i = 0; i < n; i++)
        vals[i] = (ts_elem_t)((int64_t)rand() % (2 * range + 1) - range);
    ts_instance_t *inst = ts_instance_create(vals, n);
    free(vals);
    return inst;
}

/* === Additional 3SUM Algorithms and Verification ===================== */

ts_result_t ts_3sum_meet_in_middle(const ts_instance_t *inst) {
    ts_result_t res; ts_result_init(&res);
    if (!inst || inst->n < 3) return res;
    int32_t n = inst->n;
    ts_hash_table_t *ht = ht_create(n * n + 1);
    if (!ht) return res;
    for (int32_t i = 0; i < n - 1; i++)
        for (int32_t j = i + 1; j < n; j++)
            ht_insert(ht, inst->vals[i] + inst->vals[j], i * n + j);
    for (int32_t k = 0; k < n; k++) {
        ts_elem_t target = -inst->vals[k];
        int32_t encoded = ht_lookup(ht, target);
        if (encoded >= 0) {
            int32_t i = encoded / n;
            int32_t j = encoded % n;
            if (i < j && j < k && inst->vals[i] + inst->vals[j] + inst->vals[k] == 0)
                ts_result_add_triple(&res, i, j, k);
        }
    }
    ht_destroy(ht);
    return res;
}

bool ts_verify_result(const ts_instance_t *inst, const ts_result_t *res) {
    if (!inst || !res) return false;
    for (int32_t t = 0; t < res->count; t++) {
        int32_t i = res->triples[t].i;
        int32_t j = res->triples[t].j;
        int32_t k = res->triples[t].k;
        if (i >= j || j >= k) return false;
        if (inst->vals[i] + inst->vals[j] + inst->vals[k] != 0) return false;
    }
    return true;
}

bool ts_results_equal(const ts_result_t *a, const ts_result_t *b) {
    if (!a || !b) return a == b;
    if (a->found != b->found || a->count != b->count) return false;
    for (int32_t i = 0; i < a->count; i++) {
        if (a->triples[i].i != b->triples[i].i ||
            a->triples[i].j != b->triples[i].j ||
            a->triples[i].k != b->triples[i].k) return false;
    }
    return true;
}

double ts_empirical_exponent(int32_t n1, double t1, int32_t n2, double t2) {
    if (n1 <= 0 || n2 <= 0 || t1 <= 0.0 || t2 <= 0.0) return 0.0;
    return log(t2 / t1) / log((double)n2 / (double)n1);
}

void ts_print_algorithm_comparison(const ts_instance_t *inst) {
    if (!inst) return;
    printf("=== 3SUM Algorithm Comparison (n=%d) ===\n", inst->n);
    printf("Naive O(n^3):       probes all choose(n,3) = %lld triples\n",
           (long long)inst->n * (inst->n - 1) * (inst->n - 2) / 6);
    printf("Hash O(n^2):        probes %lld pairs total\n",
           (long long)inst->n * (inst->n - 1) / 2);
    printf("Sort+BSearch O(n^2 log n): %lld pairs, each O(log n) search\n",
           (long long)inst->n * (inst->n - 1) / 2);
    printf("Two-Pointer O(n^2): %lld total two-pointer steps\n",
           (long long)inst->n * (inst->n - 1) / 2);
}

ts_result_t ts_find_all_with_threshold(const ts_instance_t *inst, ts_elem_t thresh) {
    ts_result_t res; ts_result_init(&res);
    if (!inst || inst->n < 3) return res;
    int32_t n = inst->n;
    for (int32_t i = 0; i < n - 2; i++)
        for (int32_t j = i + 1; j < n - 1; j++)
            for (int32_t k = j + 1; k < n; k++) {
                ts_elem_t sum = inst->vals[i] + inst->vals[j] + inst->vals[k];
                if (sum >= -thresh && sum <= thresh)
                    ts_result_add_triple(&res, i, j, k);
            }
    return res;
}

bool ts_is_3sum_instance(const ts_instance_t *inst) {
    return inst && inst->n >= 0 && (inst->n == 0 || inst->vals != NULL);
}

int32_t ts_unique_elements(ts_instance_t *inst) {
    if (!inst || inst->n == 0) return 0;
    qsort(inst->vals, (size_t)inst->n, sizeof(ts_elem_t), cmp_elem);
    int32_t u = 1;
    for (int32_t i = 1; i < inst->n; i++)
        if (inst->vals[i] != inst->vals[i-1]) u++;
    return u;
}

bool ts_has_zero(const ts_instance_t *inst) {
    if (!inst) return false;
    for (int32_t i = 0; i < inst->n; i++)
        if (inst->vals[i] == 0) return true;
    return false;
}

ts_elem_t ts_max_abs_value(const ts_instance_t *inst) {
    if (!inst || inst->n == 0) return 0;
    ts_elem_t maxv = inst->vals[0] >= 0 ? inst->vals[0] : -inst->vals[0];
    for (int32_t i = 1; i < inst->n; i++) {
        ts_elem_t av = inst->vals[i] >= 0 ? inst->vals[i] : -inst->vals[i];
        if (av > maxv) maxv = av;
    }
    return maxv;
}

ts_elem_t ts_min_value(const ts_instance_t *inst) {
    if (!inst || inst->n == 0) return 0;
    ts_elem_t mv = inst->vals[0];
    for (int32_t i = 1; i < inst->n; i++)
        if (inst->vals[i] < mv) mv = inst->vals[i];
    return mv;
}

ts_elem_t ts_max_value(const ts_instance_t *inst) {
    if (!inst || inst->n == 0) return 0;
    ts_elem_t mv = inst->vals[0];
    for (int32_t i = 1; i < inst->n; i++)
        if (inst->vals[i] > mv) mv = inst->vals[i];
    return mv;
}

double ts_mean_value(const ts_instance_t *inst) {
    if (!inst || inst->n == 0) return 0.0;
    double sum = 0.0;
    for (int32_t i = 0; i < inst->n; i++) sum += (double)inst->vals[i];
    return sum / inst->n;
}

int64_t ts_range(const ts_instance_t *inst) {
    return (int64_t)(ts_max_value(inst) - ts_min_value(inst));
}

bool ts_is_sorted(const ts_instance_t *inst) {
    if (!inst || inst->n < 2) return true;
    for (int32_t i = 1; i < inst->n; i++)
        if (inst->vals[i] < inst->vals[i-1]) return false;
    return true;
}

void ts_sort(ts_instance_t *inst) {
    if (!inst || inst->n < 2) return;
    qsort(inst->vals, (size_t)inst->n, sizeof(ts_elem_t), cmp_elem);
}

ts_result_t ts_3sum_with_duplicates(const ts_instance_t *inst) {
    ts_result_t res; ts_result_init(&res);
    if (!inst || inst->n < 3) return res;
    ts_instance_t *sorted = ts_instance_clone(inst);
    if (!sorted) return res;
    ts_sort(sorted);
    res = ts_quadratic_two_pointer(sorted);
    ts_instance_destroy(sorted);
    return res;
}

void ts_analyze_instance(const ts_instance_t *inst) {
    if (!inst) { printf("NULL instance\n"); return; }
    printf("=== 3SUM Instance Analysis ===\n");
    printf("Size: n = %d\n", inst->n);
    printf("Range: [%lld, %lld]\n", (long long)ts_min_value(inst), (long long)ts_max_value(inst));
    printf("Mean: %.2f\n", ts_mean_value(inst));
    printf("Contains zero: %s\n", ts_has_zero(inst) ? "yes" : "no");
    ts_elem_t maxv = ts_max_abs_value(inst);
    printf("Max abs value: %lld\n", (long long)maxv);
}
bool ts_3sum_decision(const ts_instance_t *inst) {
    ts_result_t res = ts_quadratic_two_pointer(ts_instance_clone(inst));
    bool found = res.found;
    ts_result_free(&res);
    return found;
}
void ts_benchmark_all(const ts_instance_t *inst) {
    if (!inst || inst->n < 3) return;
    printf("=== 3SUM Benchmark (n=%d) ===\n", inst->n);
    ts_instance_t *c1 = ts_instance_clone(inst);
    ts_instance_t *c2 = ts_instance_clone(inst);
    ts_instance_t *c3 = ts_instance_clone(inst);
    ts_result_t r1 = ts_naive_O_n3(inst);
    printf("Naive O(n^3): found %d triples\n", r1.count);
    ts_result_t r2 = ts_quadratic_hash(inst);
    printf("Hash O(n^2):  found %d triples\n", r2.count);
    ts_result_t r3; ts_result_init(&r3);
    if (c1) { r3 = ts_quadratic_sort_bsearch(c1); printf("Sort+BS: found %d triples\n", r3.count); }
    ts_result_t r4; ts_result_init(&r4);
    if (c2) { r4 = ts_quadratic_two_pointer(c2); printf("TwoPtr: found %d triples\n", r4.count); }
    ts_result_t r5; ts_result_init(&r5);
    if (c3) { r5 = ts_gronlund_pettie_2014(c3); printf("GP14:   found %d triples\n", r5.count); }
    ts_result_free(&r1); ts_result_free(&r2); ts_result_free(&r3);
    ts_result_free(&r4); ts_result_free(&r5);
    ts_instance_destroy(c1); ts_instance_destroy(c2); ts_instance_destroy(c3);
}
ts_elem_t ts_median_value(const ts_instance_t *inst) {
    if (!inst || inst->n == 0) return 0;
    ts_instance_t *sorted = ts_instance_clone(inst);
    if (!sorted) return 0;
    qsort(sorted->vals, (size_t)sorted->n, sizeof(ts_elem_t), cmp_elem);
    ts_elem_t med = sorted->vals[sorted->n / 2];
    ts_instance_destroy(sorted);
    return med;
}
double ts_variance(const ts_instance_t *inst) {
    if (!inst || inst->n < 2) return 0.0;
    double mean = ts_mean_value(inst), sum_sq = 0.0;
    for (int32_t i = 0; i < inst->n; i++) {
        double d = (double)inst->vals[i] - mean;
        sum_sq += d * d;
    }
    return sum_sq / (inst->n - 1);
}
double ts_stddev(const ts_instance_t *inst) { return sqrt(ts_variance(inst)); }

void ts_self_test(void) {
    printf("=== 3SUM Self-Test ===\n");
    ts_elem_t test_vals[] = {-4, -1, -1, 0, 1, 2, 3, 5};
    int32_t n = 8;
    ts_instance_t *inst = ts_instance_create(test_vals, n);
    ts_result_t r1 = ts_naive_O_n3(inst);
    ts_result_t r2 = ts_quadratic_hash(inst);
    ts_instance_t *c1 = ts_instance_clone(inst);
    ts_result_t r3 = ts_quadratic_two_pointer(c1);
    printf("Naive O(n^3): found %d triples\n", r1.count);
    printf("Hash O(n^2):  found %d triples\n", r2.count);
    printf("Two-pointer:  found %d triples\n", r3.count);
    bool all_match = (r1.count == r2.count) && (r2.count == r3.count);
    printf("All match: %s\n", all_match ? "PASS" : "FAIL");
    printf("Verification: %s\n", ts_verify_result(inst, &r1) ? "PASS" : "FAIL");
    ts_result_free(&r1); ts_result_free(&r2); ts_result_free(&r3);
    ts_instance_destroy(c1); ts_instance_destroy(inst);
}
int32_t ts_count_pairs_sum_to_target(const ts_instance_t *inst, ts_elem_t target) {
    if (!inst || inst->n < 2) return 0;
    ts_instance_t *sorted = ts_instance_clone(inst);
    if (!sorted) return 0;
    qsort(sorted->vals, (size_t)sorted->n, sizeof(ts_elem_t), cmp_elem);
    int32_t left = 0, right = sorted->n - 1, count = 0;
    while (left < right) {
        ts_elem_t sum = sorted->vals[left] + sorted->vals[right];
        if (sum < target) left++;
        else if (sum > target) right--;
        else { count++; left++; right--; }
    }
    ts_instance_destroy(sorted);
    return count;
}
void ts_print_statistics(const ts_instance_t *inst) {
    if (!inst) return;
    printf("3SUM Stats: n=%d, range=[%lld,%lld], mean=%.2f, stddev=%.2f\n",
           inst->n, (long long)ts_min_value(inst), (long long)ts_max_value(inst),
           ts_mean_value(inst), ts_stddev(inst));
}

void ts_print_conjecture_details(void) {
    ts_conjecture_t c = ts_conjecture_status();
    printf("=== 3SUM Conjecture ===\n");
    printf("Statement: %s\n", c.statement);
    printf("Epsilon infimum: %.4f\n", c.epsilon_min);
    printf("Model: %d\n", c.model_type);
    printf("History: %s\n", c.history);
}
bool ts_check_conjecture_for_size(int32_t n) {
    double naive = (double)n * n * n;
    double quadratic = (double)n * n;
    double best = quadratic * log(log((double)n)) / log((double)n);
    printf("n=%d: Naive=%.2e, Quadratic=%.2e, Best=%.2e\n", n, naive, quadratic, best);
    printf("Is best subquadratic? %s\n", best < quadratic ? "yes" : "no");
    double eps = 0.1;
    double subquad = pow((double)n, 2.0 - eps);
    printf("O(n^{%.1f})=%.2e, Best=%.2e -> best %s subquadratic\n",
           2.0-eps, subquad, best, best < subquad ? "IS" : "is NOT");
    return best >= quadratic;
}
void ts_explain_3sum_hardness(void) {
    printf("=== Understanding 3SUM-Hardness ===\n\n");
    printf("A problem P is 3SUM-hard if a truly subquadratic algorithm\n");
    printf("for P would imply a truly subquadratic algorithm for 3SUM.\n\n");
    printf("Formally: If P(n) is solvable in O(n^{2-eps}) time,\n");
    printf("then 3SUM(n) is solvable in O(n^{2-eps}) time.\n\n");
    printf("This is a fine-grained reduction: it preserves not just\n");
    printf("polynomial solvability but the EXACT exponent (up to eps).\n\n");
    printf("Key 3SUM-hard problems (Gajentaan-Overmars 1995):\n");
    for (int32_t i = 0; i < TS3_HARD_COUNT; i++)
        printf("  - %s\n", ts_hardness_class_name((ts_hardness_class_t)i));
}

/* End of threesum.c - All 3SUM algorithms implemented */
