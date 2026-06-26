/* threesum_apsp.c -- 3SUM and APSP Implementation */
#include "threesum_apsp.h"
#include "condlb.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

/* ============================================================================
 * 3SUM Implementation
 * ============================================================================ */

ThreeSumInstance* ts3_create(int capacity) {
    assert(capacity > 0);
    ThreeSumInstance* inst = (ThreeSumInstance*)malloc(sizeof(ThreeSumInstance));
    if (!inst) return NULL;
    inst->numbers = (int64_t*)malloc((size_t)capacity * sizeof(int64_t));
    if (!inst->numbers) { free(inst); return NULL; }
    inst->n = 0; inst->capacity = capacity; inst->sorted = 0;
    return inst;
}

void ts3_free(ThreeSumInstance* inst) { if(inst){free(inst->numbers);free(inst);} }

void ts3_add(ThreeSumInstance* inst, int64_t x) {
    assert(inst);
    if (inst->n >= inst->capacity) {
        int nc = inst->capacity * 2;
        int64_t* nn = (int64_t*)realloc(inst->numbers, (size_t)nc * sizeof(int64_t));
        if (!nn) return;
        inst->numbers = nn; inst->capacity = nc;
    }
    inst->numbers[inst->n++] = x;
    inst->sorted = 0;
}

static int cmp_int64(const void* a, const void* b) {
    int64_t da = *(const int64_t*)a, db = *(const int64_t*)b;
    return (da > db) - (da < db);
}

void ts3_sort(ThreeSumInstance* inst) {
    assert(inst);
    if (!inst->sorted) {
        qsort(inst->numbers, (size_t)inst->n, sizeof(int64_t), cmp_int64);
        inst->sorted = 1;
    }
}

ThreeSumSolution ts3_brute_force(const ThreeSumInstance* inst) {
    assert(inst);
    ThreeSumSolution sol = {0,0,0,0,0};
    int n = inst->n;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            for (int k = j + 1; k < n; k++)
                if (inst->numbers[i] + inst->numbers[j] + inst->numbers[k] == 0) {
                    sol.i = i; sol.j = j; sol.k = k; sol.sum = 0; sol.found = 1; return sol;
                }
    return sol;
}

ThreeSumSolution ts3_quadratic(const ThreeSumInstance* inst) {
    assert(inst);
    ThreeSumSolution sol = {0,0,0,0,0};
    int n = inst->n;
    ts3_sort((ThreeSumInstance*)inst);
    for (int i = 0; i < n - 2; i++) {
        int64_t a = inst->numbers[i];
        int left = i + 1, right = n - 1;
        while (left < right) {
            int64_t sum = a + inst->numbers[left] + inst->numbers[right];
            if (sum == 0) {
                sol.i = i; sol.j = left; sol.k = right; sol.sum = 0; sol.found = 1; return sol;
            } else if (sum < 0) left++;
            else right--;
        }
    }
    return sol;
}

ThreeSumSolution ts3_fft_small_range(const ThreeSumInstance* inst, int max_val) {
    assert(inst);
    ThreeSumSolution sol = {0,0,0,0,0};
    if (max_val <= 0) max_val = 1000;
    int range = 2 * max_val + 1;
    int* poly = (int*)calloc((size_t)range, sizeof(int));
    if (!poly) return ts3_quadratic(inst);
    for (int i = 0; i < inst->n; i++) {
        int idx = (int)(inst->numbers[i] + max_val);
        if (idx >= 0 && idx < range) poly[idx] = 1;
    }
    for (int i = 0; i < inst->n; i++) {
        int64_t a = inst->numbers[i];
        for (int j = i + 1; j < inst->n; j++) {
            int64_t need = -(a + inst->numbers[j]);
            int idx = (int)(need + max_val);
            if (idx >= 0 && idx < range && poly[idx]) {
                for (int k = j + 1; k < inst->n; k++) {
                    if (inst->numbers[k] == need) {
                        sol.i = i; sol.j = j; sol.k = k; sol.sum = 0; sol.found = 1;
                        free(poly); return sol;
                    }
                }
            }
        }
    }
    free(poly);
    return sol;
}

ThreeSumSolution ts3_randomized_hashing(const ThreeSumInstance* inst) {
    assert(inst);
    int n = inst->n;
    if (n < 3) { ThreeSumSolution sol = {0,0,0,0,0}; return sol; }
    int table_size = n * 2;
    int64_t** table = (int64_t**)calloc((size_t)table_size, sizeof(int64_t*));
    int* table_counts = (int*)calloc((size_t)table_size, sizeof(int));
    int* table_caps = (int*)calloc((size_t)table_size, sizeof(int));
    if (!table || !table_counts || !table_caps) {
        free(table); free(table_counts); free(table_caps);
        return ts3_quadratic(inst);
    }
    int found = 0;
    ThreeSumSolution sol = {0,0,0,0,0};
    srand(12345);
    uint64_t big_prime = 1000000007ULL;
    for (int i = 0; i < n && !found; i++) {
        uint64_t h = (uint64_t)(inst->numbers[i] + (1LL << 30));
        h = (h * 2654435761ULL + rand()) % big_prime % (uint64_t)table_size;
        int bucket = (int)h;
        for (int j = 0; j < table_counts[bucket] && !found; j++) {
            int idx = (int)(table[bucket][j] >> 32);
            int64_t a_plus_b = (int64_t)(table[bucket][j] & 0xFFFFFFFF);
            if (inst->numbers[i] + a_plus_b == 0) {
                for (int k = 0; k < n && !found; k++) {
                    if (k != i && k != idx && inst->numbers[idx] + inst->numbers[k] == a_plus_b) {
                        sol.i = idx; sol.j = k; sol.k = i; sol.sum = 0; sol.found = 1; found = 1;
                    }
                }
            }
        }
        if (!found) {
            for (int j = 0; j < i && !found; j++) {
                int64_t sum = inst->numbers[i] + inst->numbers[j];
                uint64_t h2 = (uint64_t)(sum + (1LL << 31));
                h2 = (h2 * 2654435761ULL + rand()) % big_prime % (uint64_t)table_size;
                int b2 = (int)h2;
                if (table_counts[b2] >= table_caps[b2]) {
                    int nc = table_caps[b2] == 0 ? 4 : table_caps[b2] * 2;
                    int64_t* nl = (int64_t*)realloc(table[b2], (size_t)nc * sizeof(int64_t));
                    if (!nl) continue;
                    table[b2] = nl; table_caps[b2] = nc;
                }
                table[b2][table_counts[b2]++] = ((int64_t)i << 32) | (sum & 0xFFFFFFFF);
            }
        }
    }
    for (int i = 0; i < table_size; i++) free(table[i]);
    free(table); free(table_counts); free(table_caps);
    return sol;
}

int ts3_exists(const ThreeSumInstance* inst) {
    ThreeSumSolution sol = ts3_quadratic(inst);
    return sol.found;
}

int ts3_count_all(const ThreeSumInstance* inst) {
    assert(inst);
    int n = inst->n, count = 0;
    ts3_sort((ThreeSumInstance*)inst);
    for (int i = 0; i < n - 2; i++) {
        int left = i + 1, right = n - 1;
        while (left < right) {
            int64_t sum = inst->numbers[i] + inst->numbers[left] + inst->numbers[right];
            if (sum == 0) {
                count++;
                int lv = left, rv = right;
                while (left + 1 < right && inst->numbers[left + 1] == inst->numbers[left]) left++;
                while (right - 1 > left && inst->numbers[right - 1] == inst->numbers[right]) right--;
                left++; right--;
            } else if (sum < 0) left++;
            else right--;
        }
    }
    return count;
}

/* ============================================================================
 * APSP Implementation
 * ============================================================================ */

ApspInstance* apsp_create(int n) {
    assert(n > 0);
    ApspInstance* inst = (ApspInstance*)malloc(sizeof(ApspInstance));
    if (!inst) return NULL;
    inst->n_vertices = n;
    inst->adj_matrix = (int64_t*)malloc((size_t)(n * n) * sizeof(int64_t));
    inst->dist = (int64_t**)malloc((size_t)n * sizeof(int64_t*));
    inst->next = (int**)malloc((size_t)n * sizeof(int*));
    if (!inst->adj_matrix || !inst->dist || !inst->next) {
        free(inst->adj_matrix); free(inst->dist); free(inst->next); free(inst); return NULL;
    }
    for (int i = 0; i < n; i++) {
        inst->dist[i] = (int64_t*)malloc((size_t)n * sizeof(int64_t));
        inst->next[i] = (int*)malloc((size_t)n * sizeof(int));
        if (!inst->dist[i] || !inst->next[i]) {
            for (int k = 0; k <= i; k++) { free(inst->dist[k]); free(inst->next[k]); }
            free(inst->adj_matrix); free(inst->dist); free(inst->next); free(inst); return NULL;
        }
    }
    int64_t INF = INT64_MAX / 2;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            inst->adj_matrix[i * n + j] = (i == j) ? 0 : INF;
            inst->dist[i][j] = (i == j) ? 0 : INF;
            inst->next[i][j] = -1;
        }
    }
    return inst;
}

void apsp_free(ApspInstance* inst) {
    if (!inst) return;
    for (int i = 0; i < inst->n_vertices; i++) { free(inst->dist[i]); free(inst->next[i]); }
    free(inst->dist); free(inst->next); free(inst->adj_matrix);
    free(inst);
}

void apsp_add_edge(ApspInstance* inst, int u, int v, int64_t w) {
    assert(inst); assert(u >= 0 && u < inst->n_vertices); assert(v >= 0 && v < inst->n_vertices);
    int n = inst->n_vertices;
    inst->adj_matrix[u * n + v] = w;
    inst->dist[u][v] = w;
    inst->next[u][v] = v;
}

void apsp_floyd_warshall(ApspInstance* inst) {
    assert(inst);
    int n = inst->n_vertices;
    int64_t INF = INT64_MAX / 2;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            inst->dist[i][j] = inst->adj_matrix[i * n + j];
            inst->next[i][j] = (inst->adj_matrix[i * n + j] < INF && i != j) ? j : -1;
        }
    for (int k = 0; k < n; k++)
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                if (inst->dist[i][k] < INF && inst->dist[k][j] < INF &&
                    inst->dist[i][k] + inst->dist[k][j] < inst->dist[i][j]) {
                    inst->dist[i][j] = inst->dist[i][k] + inst->dist[k][j];
                    inst->next[i][j] = inst->next[i][k];
                }
}

void apsp_dijkstra_all(ApspInstance* inst) {
    assert(inst);
    int n = inst->n_vertices;
    int64_t INF = INT64_MAX / 2;
    for (int src = 0; src < n; src++) {
        int* visited = (int*)calloc((size_t)n, sizeof(int));
        if (!visited) continue;
        for (int i = 0; i < n; i++) { inst->dist[src][i] = INF; inst->next[src][i] = -1; }
        inst->dist[src][src] = 0;
        for (int iter = 0; iter < n; iter++) {
            int u = -1; int64_t min_d = INF;
            for (int v = 0; v < n; v++)
                if (!visited[v] && inst->dist[src][v] < min_d) { min_d = inst->dist[src][v]; u = v; }
            if (u == -1) break;
            visited[u] = 1;
            for (int v = 0; v < n; v++) {
                int64_t w = inst->adj_matrix[u * n + v];
                if (w < INF && inst->dist[src][u] + w < inst->dist[src][v]) {
                    inst->dist[src][v] = inst->dist[src][u] + w;
                    inst->next[src][v] = (inst->next[src][u] == -1) ? v : inst->next[src][u];
                }
            }
        }
        free(visited);
    }
}

void apsp_min_plus_product(ApspInstance* inst) {
    apsp_floyd_warshall(inst);
}

void apsp_small_weights(ApspInstance* inst, int max_weight) {
    assert(inst && max_weight > 0);
    apsp_floyd_warshall(inst);
}

int apsp_has_negative_cycle(const ApspInstance* inst) {
    assert(inst);
    int n = inst->n_vertices;
    int64_t INF = INT64_MAX / 2;
    for (int i = 0; i < n; i++)
        if (inst->dist[i][i] < 0) return 1;
    return 0;
}

int apsp_reconstruct_path(const ApspInstance* inst, int u, int v, int* path, int max_len) {
    assert(inst && path && max_len > 0);
    if (inst->dist[u][v] >= INT64_MAX / 4) return -1;
    int len = 0;
    path[len++] = u;
    while (u != v && len < max_len) {
        u = inst->next[u][v];
        if (u == -1) return -1;
        path[len++] = u;
    }
    return len;
}

/* ============================================================================
 * CLB Theorems and Applications
 * ============================================================================ */

int ts3_check_conjecture(double algo_exponent) {
    return (algo_exponent < 1.999) ? 1 : 0;
}

double apsp_subcubic_equivalence_exponent(void) { return 3.0; }

int apsp_subcubic_class_size(void) { return 32; }

int geom_three_collinear(const int64_t* x, const int64_t* y, int n) {
    assert(x && y); assert(n >= 0);
    if (n < 3) return 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            for (int k = j + 1; k < n; k++) {
                int64_t cross = (x[j] - x[i]) * (y[k] - y[i]) - (y[j] - y[i]) * (x[k] - x[i]);
                if (cross == 0) return 1;
            }
    return 0;
}

int64_t apsp_graph_radius(const ApspInstance* inst) {
    assert(inst);
    int n = inst->n_vertices;
    int64_t INF = INT64_MAX / 2, radius = INF;
    for (int v = 0; v < n; v++) {
        int64_t eccentricity = 0;
        for (int u = 0; u < n; u++) {
            if (inst->dist[v][u] > eccentricity) eccentricity = inst->dist[v][u];
        }
        if (eccentricity < radius) radius = eccentricity;
    }
    return (radius >= INF) ? -1 : radius;
}

void apsp_betweenness_centrality(const ApspInstance* inst, double* centrality) {
    assert(inst && centrality);
    int n = inst->n_vertices;
    for (int v = 0; v < n; v++) centrality[v] = 0.0;
    for (int s = 0; s < n; s++) {
        for (int t = 0; t < n; t++) {
            if (s == t) continue;
            for (int v = 0; v < n; v++) {
                if (v == s || v == t) continue;
                int64_t d_st = inst->dist[s][t];
                int64_t d_sv = inst->dist[s][v];
                int64_t d_vt = inst->dist[v][t];
                if (d_st < INT64_MAX/4 && d_sv < INT64_MAX/4 && d_vt < INT64_MAX/4 &&
                    d_sv + d_vt == d_st) {
                    centrality[v] += 1.0;
                }
            }
        }
    }
    for (int v = 0; v < n; v++) centrality[v] /= 2.0;
}

int64_t geom_min_area_triangle(const int64_t* x, const int64_t* y, int n) {
    assert(x && y);
    if (n < 3) return INT64_MAX;
    int64_t min_area = INT64_MAX;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            for (int k = j + 1; k < n; k++) {
                int64_t area = (x[j] - x[i]) * (y[k] - y[i]) - (y[j] - y[i]) * (x[k] - x[i]);
                if (area < 0) area = -area;
                if (area > 0 && area < min_area) min_area = area;
            }
    return (min_area == INT64_MAX) ? 0 : min_area;
}
