/* seth_bound.c -- SETH-based Conditional Lower Bounds: Implementation */
#include "seth_bound.h"
#include "condlb.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

/* ============================================================================
 * L1: SAT Instance Management
 * ============================================================================ */

SatInstance* sat_create(int n_vars, int n_clauses, int k) {
    assert(n_vars > 0 && n_clauses > 0 && k >= 1);
    SatInstance* inst = (SatInstance*)malloc(sizeof(SatInstance));
    if (!inst) return NULL;
    inst->n_vars = n_vars;
    inst->n_clauses = 0;
    inst->max_clause_len = k;
    inst->capacity = n_clauses;
    inst->clauses = (int**)calloc((size_t)n_clauses, sizeof(int*));
    inst->clause_sizes = (int*)calloc((size_t)n_clauses, sizeof(int));
    if (!inst->clauses || !inst->clause_sizes) {
        free(inst->clauses); free(inst->clause_sizes); free(inst); return NULL;
    }
    return inst;
}

void sat_free(SatInstance* inst) {
    if (!inst) return;
    if (inst->clauses) {
        for (int i = 0; i < inst->n_clauses; i++) free(inst->clauses[i]);
        free(inst->clauses);
    }
    free(inst->clause_sizes);
    free(inst);
}

int sat_add_clause(SatInstance* inst, const int* literals, int len) {
    assert(inst != NULL && literals != NULL && len >= 1);
    if (inst->n_clauses >= inst->capacity) return -1;
    int idx = inst->n_clauses++;
    inst->clause_sizes[idx] = len;
    inst->clauses[idx] = (int*)malloc((size_t)len * sizeof(int));
    if (!inst->clauses[idx]) return -1;
    memcpy(inst->clauses[idx], literals, (size_t)len * sizeof(int));
    if (len > inst->max_clause_len) inst->max_clause_len = len;
    return idx;
}

int sat_check(const SatInstance* inst, const int* assignment) {
    assert(inst != NULL && assignment != NULL);
    for (int i = 0; i < inst->n_clauses; i++) {
        int clause_sat = 0;
        for (int j = 0; j < inst->clause_sizes[i]; j++) {
            int lit = inst->clauses[i][j];
            int var_idx = (lit > 0) ? lit - 1 : -lit - 1;
            int val = (lit > 0) ? 1 : 0;
            if (assignment[var_idx] == val) { clause_sat = 1; break; }
        }
        if (!clause_sat) return 0;
    }
    return 1;
}

/* ============================================================================
 * L5: SAT Solving Algorithms
 * ============================================================================ */

int sat_brute_force(const SatInstance* inst, int* assignment) {
    assert(inst != NULL);
    int n = inst->n_vars;
    int total = 1 << n;
    int* try_assign = assignment ? assignment : (int*)malloc((size_t)n * sizeof(int));
    int need_free = (assignment == NULL);
    if (!try_assign) return -1;
    for (int mask = 0; mask < total; mask++) {
        for (int v = 0; v < n; v++) try_assign[v] = (mask >> v) & 1;
        if (sat_check(inst, try_assign)) {
            if (need_free) free(try_assign);
            return 1;
        }
    }
    if (need_free) free(try_assign);
    return 0;
}

static int dpll_recursive(SatInstance* inst, int* assign, int depth) {
    int n = inst->n_vars;
    int all_sat = 1;
    for (int i = 0; i < inst->n_clauses; i++) {
        int clause_sat = 0, unassigned_count = 0, unassigned_lit = 0;
        for (int j = 0; j < inst->clause_sizes[i]; j++) {
            int lit = inst->clauses[i][j];
            int var_idx = (lit > 0) ? lit - 1 : -lit - 1;
            int val = (lit > 0) ? 1 : 0;
            if (assign[var_idx] == -1) { unassigned_count++; unassigned_lit = lit; }
            else if (assign[var_idx] == val) { clause_sat = 1; break; }
        }
        if (!clause_sat) {
            if (unassigned_count == 0) return 0;
            if (unassigned_count == 1) {
                int vi = (unassigned_lit > 0) ? unassigned_lit - 1 : -unassigned_lit - 1;
                int val = (unassigned_lit > 0) ? 1 : 0;
                if (assign[vi] == 1 - val) return 0;
                assign[vi] = val;
                int result = dpll_recursive(inst, assign, depth + 1);
                if (result) return result;
                assign[vi] = -1;
            }
            all_sat = 0;
        }
    }
    if (all_sat) return 1;
    for (int v = 0; v < n; v++) {
        if (assign[v] != -1) continue;
        int pos_count = 0, neg_count = 0;
        for (int i = 0; i < inst->n_clauses; i++)
            for (int j = 0; j < inst->clause_sizes[i]; j++) {
                int lit = inst->clauses[i][j];
                int vi = (lit > 0) ? lit - 1 : -lit - 1;
                if (vi == v) { if (lit > 0) pos_count++; else neg_count++; }
            }
        if (pos_count > 0 && neg_count == 0) {
            assign[v] = 1; int r = dpll_recursive(inst, assign, depth + 1);
            assign[v] = -1; return r;
        }
        if (neg_count > 0 && pos_count == 0) {
            assign[v] = 0; int r = dpll_recursive(inst, assign, depth + 1);
            assign[v] = -1; return r;
        }
    }
    int best_var = -1, best_count = -1;
    for (int v = 0; v < n; v++) {
        if (assign[v] != -1) continue;
        int count = 0;
        for (int i = 0; i < inst->n_clauses; i++)
            for (int j = 0; j < inst->clause_sizes[i]; j++) {
                int lit = inst->clauses[i][j];
                if (((lit > 0) ? lit - 1 : -lit - 1) == v) count++;
            }
        if (count > best_count) { best_count = count; best_var = v; }
    }
    if (best_var == -1) return 0;
    assign[best_var] = 1;
    if (dpll_recursive(inst, assign, depth + 1)) return 1;
    assign[best_var] = 0;
    if (dpll_recursive(inst, assign, depth + 1)) return 1;
    assign[best_var] = -1;
    return 0;
}

int sat_dpll(const SatInstance* inst, int* assignment) {
    assert(inst != NULL);
    int n = inst->n_vars;
    int* assign = assignment ? assignment : (int*)malloc((size_t)n * sizeof(int));
    int need_free = (assignment == NULL);
    if (!assign) return -1;
    for (int i = 0; i < n; i++) assign[i] = -1;
    int result = dpll_recursive((SatInstance*)inst, assign, 0);
    if (!result && need_free) free(assign);
    return result;
}

int sat_schoening(const SatInstance* inst, int* assignment, int max_tries) {
    assert(inst != NULL);
    int n = inst->n_vars;
    if (max_tries <= 0) {
        double k = (double)inst->max_clause_len;
        double base = 2.0 * (k - 1.0) / k;
        max_tries = (int)(pow(base, n) * 10.0);
        if (max_tries < 100) max_tries = 100;
    }
    int* assign = assignment ? assignment : (int*)malloc((size_t)n * sizeof(int));
    int need_free = (assignment == NULL);
    if (!assign) return -1;
    srand((unsigned)time(NULL));
    for (int t = 0; t < max_tries; t++) {
        for (int v = 0; v < n; v++) assign[v] = rand() % 2;
        for (int step = 0; step < 3 * n; step++) {
            int unsat_clause = -1;
            for (int i = 0; i < inst->n_clauses; i++) {
                int sat = 0;
                for (int j = 0; j < inst->clause_sizes[i]; j++) {
                    int lit = inst->clauses[i][j];
                    int vi = (lit > 0) ? lit - 1 : -lit - 1;
                    if (assign[vi] == ((lit > 0) ? 1 : 0)) { sat = 1; break; }
                }
                if (!sat) { unsat_clause = i; break; }
            }
            if (unsat_clause < 0) { if (need_free) free(assign); return 1; }
            int rand_idx = rand() % inst->clause_sizes[unsat_clause];
            int lit = inst->clauses[unsat_clause][rand_idx];
            int vi = (lit > 0) ? lit - 1 : -lit - 1;
            assign[vi] = 1 - assign[vi];
        }
    }
    if (need_free) free(assign);
    return 0;
}

int sat_ppsz(const SatInstance* inst, int* assignment) {
    assert(inst != NULL);
    int n = inst->n_vars;
    int* assign = assignment ? assignment : (int*)malloc((size_t)n * sizeof(int));
    int need_free = (assignment == NULL);
    if (!assign) return -1;
    int trials = (int)(pow(2.0, (double)n * 0.39) + 100);
    if (trials > 500000) trials = 500000;
    srand((unsigned)time(NULL));
    for (int trial = 0; trial < trials; trial++) {
        int* pi = (int*)malloc((size_t)n * sizeof(int));
        if (!pi) continue;
        for (int i = 0; i < n; i++) pi[i] = i;
        for (int i = n - 1; i > 0; i--) {
            int j = rand() % (i + 1); int tmp = pi[i]; pi[i] = pi[j]; pi[j] = tmp;
        }
        for (int v = 0; v < n; v++) assign[v] = -1;
        for (int idx = 0; idx < n; idx++) {
            int v = pi[idx]; assign[v] = rand() % 2;
            int changed = 1;
            while (changed) {
                changed = 0;
                for (int i = 0; i < inst->n_clauses; i++) {
                    int uc = 0, ul = 0, cs = 0;
                    for (int j = 0; j < inst->clause_sizes[i]; j++) {
                        int lit = inst->clauses[i][j];
                        int vi = (lit > 0) ? lit - 1 : -lit - 1;
                        if (assign[vi] == -1) { uc++; ul = lit; }
                        else if (assign[vi] == ((lit > 0) ? 1 : 0)) { cs = 1; break; }
                    }
                    if (!cs && uc == 1) {
                        int vi = (ul > 0) ? ul - 1 : -ul - 1;
                        assign[vi] = (ul > 0) ? 1 : 0;
                        changed = 1;
                    }
                }
            }
        }
        for (int v = 0; v < n; v++) if (assign[v] == -1) assign[v] = 0;
        if (sat_check(inst, assign)) { free(pi); if (need_free) free(assign); return 1; }
        free(pi);
    }
    if (need_free) free(assign);
    return 0;
}

/* ============================================================================
 * L4: SETH Theorems and Bounds
 * ============================================================================ */

double seth_lower_bound_for_k(int k) {
    assert(k >= 3);
    static const double sk[] = {0.0,0.0,0.0, 0.386, 0.559, 0.659, 0.726, 0.775};
    if (k <= 7) return sk[k];
    return 1.0 - log((double)k) / (double)k;
}

double seth_circuit_sat_bound(int depth_d, int n_vars) {
    assert(depth_d >= 1 && n_vars > 0);
    double ln = log((double)n_vars);
    double denom = 1.0;
    for (int i = 1; i < depth_d; i++) denom *= ln;
    return 1.0 - 1.0 / denom;
}

double seth_parameterized_bound(int n, int parameter_p) {
    assert(n > 0 && parameter_p > 0);
    if (parameter_p <= 0) return 0.0;
    return (double)parameter_p * 0.1;
}

/* ============================================================================
 * L5: SETH-to-Problem Reductions
 * ============================================================================ */

void seth_to_ov(const SatInstance* sat_inst,
                int** out_vec_a, int* n_a, int* dim_a,
                int** out_vec_b, int* n_b, int* dim_b) {
    assert(sat_inst != NULL);
    int n = sat_inst->n_vars, n_half = n / 2, m = sat_inst->n_clauses;
    *dim_a = m; *dim_b = m;
    *n_a = 1 << n_half; *n_b = 1 << (n - n_half);
    int* va = (int*)calloc((size_t)(*n_a * m), sizeof(int));
    int* vb = (int*)calloc((size_t)(*n_b * m), sizeof(int));
    assert(va && vb);
    for (int au = 0; au < *n_a; au++) {
        int* p = (int*)calloc((size_t)n, sizeof(int));
        for (int v = 0; v < n_half; v++) p[v] = (au >> v) & 1;
        for (int c = 0; c < m; c++) {
            int cf = 1;
            for (int j = 0; j < sat_inst->clause_sizes[c]; j++) {
                int lit = sat_inst->clauses[c][j];
                int vi = (lit > 0) ? lit - 1 : -lit - 1;
                if (vi < n_half) { if (p[vi] == ((lit > 0) ? 1 : 0)) { cf = 0; break; } }
                else { cf = 0; break; }
            }
            va[au * m + c] = cf ? 1 : 0;
        }
        free(p);
    }
    for (int av = 0; av < *n_b; av++) {
        int* p = (int*)calloc((size_t)n, sizeof(int));
        for (int v = n_half; v < n; v++) p[v] = (av >> (v - n_half)) & 1;
        for (int c = 0; c < m; c++) {
            int cf = 1;
            for (int j = 0; j < sat_inst->clause_sizes[c]; j++) {
                int lit = sat_inst->clauses[c][j];
                int vi = (lit > 0) ? lit - 1 : -lit - 1;
                if (vi >= n_half) { if (p[vi] == ((lit > 0) ? 1 : 0)) { cf = 0; break; } }
                else { cf = 0; break; }
            }
            vb[av * m + c] = cf ? 1 : 0;
        }
        free(p);
    }
    *out_vec_a = va; *out_vec_b = vb;
}

void seth_to_edit_distance(const SatInstance* sat_inst,
                           char** out_s1, int* len1,
                           char** out_s2, int* len2) {
    assert(sat_inst != NULL);
    int n = sat_inst->n_vars, m = sat_inst->n_clauses, nh = n / 2, N = 1 << nh;
    *len1 = N * (m + 2); *out_s1 = (char*)calloc((size_t)(*len1 + 1), 1);
    *len2 = N * (m + 2); *out_s2 = (char*)calloc((size_t)(*len2 + 1), 1);
    if (!*out_s1 || !*out_s2) { free(*out_s1); free(*out_s2); return; }
    for (int i = 0; i < N; i++) {
        int base = i * (m + 2);
        (*out_s1)[base] = '#'; (*out_s2)[base] = '#';
        for (int b = 0; b < nh; b++) {
            int bi = base + 1 + b;
            if (bi < *len1) {
                (*out_s1)[bi] = ((i >> b) & 1) ? '1' : '0';
                (*out_s2)[bi] = ((i >> b) & 1) ? '1' : '0';
            }
        }
        int* p = (int*)calloc((size_t)n, sizeof(int));
        for (int v = 0; v < nh; v++) p[v] = (i >> v) & 1;
        for (int c = 0; c < m; c++) {
            int cs = 0;
            for (int j = 0; j < sat_inst->clause_sizes[c]; j++) {
                int lit = sat_inst->clauses[c][j];
                int vi = (lit > 0) ? lit - 1 : -lit - 1;
                if (vi < nh && p[vi] == ((lit > 0) ? 1 : 0)) cs = 1;
            }
            int idx = base + 1 + nh + c;
            if (idx < *len1) { (*out_s1)[idx] = cs ? 'Y' : 'N'; (*out_s2)[idx] = cs ? 'Y' : 'N'; }
        }
        free(p);
    }
}

void seth_to_hitting_set(const SatInstance* sat_inst,
                         int** out_universe, int* universe_size,
                         int*** out_sets, int* n_sets, int** set_sizes) {
    assert(sat_inst != NULL);
    int n = sat_inst->n_vars, nh = n / 2, nh2 = n - nh, m = sat_inst->n_clauses;
    *universe_size = (1 << nh) + (1 << nh2);
    *out_universe = (int*)malloc((size_t)(*universe_size) * sizeof(int));
    if (!*out_universe) return;
    for (int i = 0; i < *universe_size; i++) (*out_universe)[i] = i;
    *n_sets = m + 2;
    *set_sizes = (int*)calloc((size_t)(*n_sets), sizeof(int));
    *out_sets = (int**)malloc((size_t)(*n_sets) * sizeof(int*));
    if (!*set_sizes || !*out_sets) { free(*out_universe); free(*set_sizes); free(*out_sets); return; }
    for (int c = 0; c < m; c++) {
        int sz = 0;
        for (int au = 0; au < (1 << nh); au++) {
            int cas = 0;
            for (int j = 0; j < sat_inst->clause_sizes[c]; j++) {
                int lit = sat_inst->clauses[c][j];
                int vi = (lit > 0) ? lit - 1 : -lit - 1;
                if (vi < nh && ((au >> vi) & 1) == ((lit > 0) ? 1 : 0)) { cas = 1; break; }
            }
            if (!cas) sz++;
        }
        (*set_sizes)[c] = sz;
        (*out_sets)[c] = (int*)malloc((size_t)sz * sizeof(int));
        int idx = 0;
        for (int au = 0; au < (1 << nh); au++) {
            int cas = 0;
            for (int j = 0; j < sat_inst->clause_sizes[c]; j++) {
                int lit = sat_inst->clauses[c][j];
                int vi = (lit > 0) ? lit - 1 : -lit - 1;
                if (vi < nh && ((au >> vi) & 1) == ((lit > 0) ? 1 : 0)) { cas = 1; break; }
            }
            if (!cas && idx < sz) (*out_sets)[c][idx++] = au;
        }
    }
    (*set_sizes)[m] = 1 << nh;
    (*out_sets)[m] = (int*)malloc((size_t)(1 << nh) * sizeof(int));
    for (int i = 0; i < (1 << nh); i++) (*out_sets)[m][i] = i;
    (*set_sizes)[m + 1] = 1 << nh2;
    (*out_sets)[m + 1] = (int*)malloc((size_t)(1 << nh2) * sizeof(int));
    for (int i = 0; i < (1 << nh2); i++) (*out_sets)[m + 1][i] = (1 << nh) + i;
}

/* ============================================================================
 * L7-L8: SETH Status and Analysis
 * ============================================================================ */

int seth_check_refutation(double algorithm_exponent,
                          const char* problem_name,
                          const char* reduction_ref) {
    assert(problem_name != NULL);
    (void)reduction_ref;
    double threshold = 2.0;
    if (strstr(problem_name, "Edit") || strstr(problem_name, "LCS") ||
        strstr(problem_name, "Frechet") || strstr(problem_name, "DTW") ||
        strstr(problem_name, "OV") || strstr(problem_name, "Diameter") ||
        strstr(problem_name, "Hitting") || strstr(problem_name, "Regex"))
        threshold = 2.0;
    return (algorithm_exponent < threshold - 1e-9) ? 1 : 0;
}

double seth_threshold(const char* problem_name) {
    assert(problem_name != NULL);
    if (strstr(problem_name, "OV") || strstr(problem_name, "Edit") ||
        strstr(problem_name, "LCS") || strstr(problem_name, "Frechet") ||
        strstr(problem_name, "DTW") || strstr(problem_name, "Diameter") ||
        strstr(problem_name, "Hitting")) return 2.0;
    return 2.0;
}

void seth_status_report(void) {
    printf("=== SETH Status Report ===\n\n");
    printf("SETH (Strong Exponential Time Hypothesis):\n");
    printf("  Proposed by: Impagliazzo and Paturi (2001)\n");
    printf("  Statement: lim_{k->inf} s_k = 1\n");
    printf("  Status: UNREFUTED (widely believed)\n\n");
    printf("Best known upper bounds for s_k:\n");
    printf("  k=3: s_3 <= %.4f  (PPSZ O(1.307^n))\n", seth_lower_bound_for_k(3));
    printf("  k=4: s_4 <= %.4f  (PPSZ O(1.559^n))\n", seth_lower_bound_for_k(4));
    printf("  k=5: s_5 <= %.4f  (PPSZ O(1.659^n))\n", seth_lower_bound_for_k(5));
    printf("  k=6: s_6 <= %.4f  (PPSZ O(1.726^n))\n", seth_lower_bound_for_k(6));
    printf("  k=7: s_7 <= %.4f  (PPSZ O(1.775^n))\n", seth_lower_bound_for_k(7));
    printf("  As k -> inf: s_k -> 1 (conjectured)\n\n");
    printf("Key SETH-based conditional lower bounds:\n");
    printf("  - Orthogonal Vectors (Williams, 2005)\n");
    printf("  - Edit Distance (Backurs-Indyk, STOC 2015)\n");
    printf("  - LCS (Abboud-Backurs-Williams, FOCS 2015)\n");
    printf("  - Frechet Distance (Bringmann, FOCS 2014)\n");
    printf("  - DTW (ABB FOCS 2015; BK FOCS 2015)\n");
    printf("  - Diameter (Roditty-Williams, ICALP 2013)\n");
    printf("  - Hitting Set (Abboud et al., STOC 2016)\n");
    printf("  - Regex Matching (Backurs-Indyk, FOCS 2016)\n");
    printf("  - Subgraph Isomorphism (BGL, SODA 2017)\n\n");
    printf("NSETH (Nondeterministic SETH):\n");
    printf("  Proposed by: Carmosino et al. (2016)\n");
    printf("  Status: UNREFUTED\n");
}
double nseth_lower_bound(const char* problem_name) {
    assert(problem_name != NULL);
    if (strstr(problem_name, "TAUT") || strstr(problem_name, "tautology"))
        return 0.999;
    return 0.0;
}

int nseth_consistency_check(double nondet_exponent, double det_exponent) {
    if (nondet_exponent < 0.999) return 0;
    if (det_exponent < nondet_exponent) return 0;
    return 1;
}
