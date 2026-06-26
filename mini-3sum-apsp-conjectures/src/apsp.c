/* apsp.c - All-Pairs Shortest Paths algorithms */
#include "apsp.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>

#define DMIN(a,b) ((a)<(b)?(a):(b))

apsp_adjacency_t *apsp_adjacency_create(int32_t n) {
    apsp_adjacency_t *adj = malloc(sizeof(apsp_adjacency_t));
    if (!adj) return NULL;
    adj->n = n;
    adj->weights = malloc((size_t)(n * n) * sizeof(double));
    if (!adj->weights) { free(adj); return NULL; }
    for (int32_t i = 0; i < n * n; i++)
        adj->weights[i] = (i / n == i % n) ? 0.0 : APSP_INF;
    return adj;
}

void apsp_adjacency_destroy(apsp_adjacency_t *adj) {
    if (!adj) return;
    free(adj->weights);
    free(adj);
}

void apsp_adjacency_set_edge(apsp_adjacency_t *adj, int32_t u, int32_t v, double w) {
    if (!adj || u < 0 || u >= adj->n || v < 0 || v >= adj->n) return;
    adj->weights[u * adj->n + v] = w;
}

apsp_graph_t *apsp_graph_create(int32_t n) {
    apsp_graph_t *g = malloc(sizeof(apsp_graph_t));
    if (!g) return NULL;
    g->n = n;
    g->has_negative_cycle = false;
    g->dist = malloc((size_t)(n * n) * sizeof(double));
    g->next = malloc((size_t)(n * n) * sizeof(int32_t));
    if (!g->dist || !g->next) {
        free(g->dist); free(g->next); free(g); return NULL;
    }
    for (int32_t i = 0; i < n; i++) {
        for (int32_t j = 0; j < n; j++) {
            g->dist[i * n + j] = (i == j) ? 0.0 : APSP_INF;
            g->next[i * n + j] = -1;
        }
    }
    return g;
}

void apsp_graph_destroy(apsp_graph_t *g) {
    if (!g) return;
    free(g->dist); free(g->next);
    free(g);
}

void apsp_print_graph(const apsp_graph_t *g) {
    if (!g) { printf("NULL graph\n"); return; }
    printf("APSP Graph: n=%d vertices\n", g->n);
    if (g->has_negative_cycle) printf("WARNING: negative cycle detected!\n");
    int32_t display = g->n < 8 ? g->n : 8;
    printf("Distance matrix (%dx%d shown):\n", display, display);
    for (int32_t i = 0; i < display; i++) {
        for (int32_t j = 0; j < display; j++) {
            double d = g->dist[i * g->n + j];
            if (d >= APSP_INF / 2.0) printf("  INF");
            else printf(" %4.1f", d);
        }
        printf("\n");
    }
    if (g->n > 8) printf("  ... (%d more rows/cols)\n", g->n - 8);
}

apsp_adjacency_t *apsp_generate_random_graph(int32_t n, double density, double max_weight, uint64_t seed) {
    if (seed == 0) seed = (uint64_t)time(NULL);
    srand((unsigned int)seed);
    apsp_adjacency_t *adj = apsp_adjacency_create(n);
    if (!adj) return NULL;
    for (int32_t i = 0; i < n; i++) {
        for (int32_t j = 0; j < n; j++) {
            if (i != j && ((double)rand() / RAND_MAX) < density) {
                double w = ((double)rand() / RAND_MAX) * max_weight + 0.1;
                adj->weights[i * n + j] = w;
            }
        }
    }
    return adj;
}

apsp_conjecture_t apsp_conjecture_status(void) {
    apsp_conjecture_t c;
    c.statement = "For every epsilon > 0, there is no algorithm that solves "
                  "APSP on n-vertex graphs with arbitrary real edge weights "
                  "in O(n^{3-epsilon}) time. (Vassilevska-Williams 2009)";
    c.epsilon_min = 0.0;
    c.history = "2009: Subcubic equivalence established; "
                "APSAP <-> Min-Plus Product <-> Negative Triangle; "
                "Seidel (1995): O(M(n) log n) for unweighted undirected; "
                "Zwick (2002): O(n^{2.575}) for directed with bounded weights; "
                "Williams (2014): O(n^3 / 2^{Omega(sqrt(log n))}); "
                "O(n^{3-epsilon}) barrier remains for general case.";
    return c;
}

apsp_graph_t *apsp_floyd_warshall(const apsp_adjacency_t *adj) {
    if (!adj) return NULL;
    int32_t n = adj->n;
    apsp_graph_t *g = apsp_graph_create(n);
    if (!g) return NULL;
    for (int32_t i = 0; i < n; i++) {
        for (int32_t j = 0; j < n; j++) {
            g->dist[i * n + j] = adj->weights[i * n + j];
            if (i != j && adj->weights[i * n + j] < APSP_INF / 2.0)
                g->next[i * n + j] = j;
        }
    }
    for (int32_t k = 0; k < n; k++) {
        for (int32_t i = 0; i < n; i++) {
            if (g->dist[i * n + k] >= APSP_INF / 2.0) continue;
            for (int32_t j = 0; j < n; j++) {
                if (g->dist[k * n + j] >= APSP_INF / 2.0) continue;
                double nd = g->dist[i * n + k] + g->dist[k * n + j];
                if (nd < g->dist[i * n + j]) {
                    g->dist[i * n + j] = nd;
                    g->next[i * n + j] = g->next[i * n + k];
                }
            }
        }
    }
    for (int32_t i = 0; i < n; i++)
        if (g->dist[i * n + i] < 0.0)
            g->has_negative_cycle = true;
    return g;
}

/* Bellman-Ford single-source shortest paths, returns true if no negative cycle */
static bool bellman_ford(const apsp_adjacency_t *adj, int32_t src, double *dist) {
    int32_t n = adj->n;
    for (int32_t i = 0; i < n; i++) dist[i] = APSP_INF;
    dist[src] = 0.0;
    for (int32_t iter = 0; iter < n - 1; iter++) {
        bool changed = false;
        for (int32_t u = 0; u < n; u++) {
            if (dist[u] >= APSP_INF / 2.0) continue;
            for (int32_t v = 0; v < n; v++) {
                double w = adj->weights[u * n + v];
                if (w >= APSP_INF / 2.0) continue;
                if (dist[u] + w < dist[v]) {
                    dist[v] = dist[u] + w;
                    changed = true;
                }
            }
        }
        if (!changed) break;
    }
    for (int32_t u = 0; u < n; u++) {
        if (dist[u] >= APSP_INF / 2.0) continue;
        for (int32_t v = 0; v < n; v++) {
            double w = adj->weights[u * n + v];
            if (w >= APSP_INF / 2.0) continue;
            if (dist[u] + w < dist[v] - 1e-12) return false;
        }
    }
    return true;
}

/* Dijkstra single-source for non-negative weights, O(n^2) */
static void dijkstra(const apsp_adjacency_t *adj, int32_t src,
                     double *dist, int32_t *prev) {
    int32_t n = adj->n;
    bool *visited = calloc((size_t)n, sizeof(bool));
    for (int32_t i = 0; i < n; i++) { dist[i] = APSP_INF; prev[i] = -1; }
    dist[src] = 0.0;
    for (int32_t iter = 0; iter < n; iter++) {
        int32_t u = -1;
        double min_d = APSP_INF;
        for (int32_t i = 0; i < n; i++) {
            if (!visited[i] && dist[i] < min_d) {
                min_d = dist[i]; u = i;
            }
        }
        if (u < 0 || min_d >= APSP_INF / 2.0) break;
        visited[u] = true;
        for (int32_t v = 0; v < n; v++) {
            double w = adj->weights[u * n + v];
            if (w >= APSP_INF / 2.0 || visited[v]) continue;
            double nd = dist[u] + w;
            if (nd < dist[v]) { dist[v] = nd; prev[v] = u; }
        }
    }
    free(visited);
}

apsp_graph_t *apsp_dijkstra_all(const apsp_adjacency_t *adj) {
    if (!adj) return NULL;
    int32_t n = adj->n;
    apsp_graph_t *g = apsp_graph_create(n);
    if (!g) return NULL;
    double *dist = malloc((size_t)n * sizeof(double));
    int32_t *prev = malloc((size_t)n * sizeof(int32_t));
    if (!dist || !prev) { free(dist); free(prev); apsp_graph_destroy(g); return NULL; }
    for (int32_t s = 0; s < n; s++) {
        dijkstra(adj, s, dist, prev);
        for (int32_t t = 0; t < n; t++) {
            g->dist[s * n + t] = dist[t];
            g->next[s * n + t] = (dist[t] < APSP_INF / 2.0 && t != s) ? prev[t] : -1;
        }
    }
    free(dist); free(prev);
    return g;
}

apsp_graph_t *apsp_johnson(apsp_adjacency_t *adj) {
    if (!adj) return NULL;
    int32_t n = adj->n;
    double *h = malloc((size_t)n * sizeof(double));
    if (!h) return NULL;
    if (!bellman_ford(adj, 0, h)) {
        bool ok = false;
        for (int32_t s = 1; s < n && !ok; s++)
            ok = bellman_ford(adj, s, h);
        if (!ok) {
            free(h);
            apsp_graph_t *g = apsp_graph_create(n);
            if (g) g->has_negative_cycle = true;
            return g;
        }
    }
    apsp_adjacency_t *rew = apsp_adjacency_create(n);
    if (!rew) { free(h); return NULL; }
    for (int32_t u = 0; u < n; u++)
        for (int32_t v = 0; v < n; v++)
            if (adj->weights[u * n + v] < APSP_INF / 2.0)
                rew->weights[u * n + v] = adj->weights[u * n + v] + h[u] - h[v];
    apsp_graph_t *g = apsp_dijkstra_all(rew);
    if (g) {
        for (int32_t u = 0; u < n; u++)
            for (int32_t v = 0; v < n; v++)
                if (g->dist[u * n + v] < APSP_INF / 2.0)
                    g->dist[u * n + v] = g->dist[u * n + v] - h[u] + h[v];
    }
    apsp_adjacency_destroy(rew); free(h);
    return g;
}

bool apsp_has_negative_cycle(const apsp_adjacency_t *adj) {
    if (!adj) return false;
    double *dist = malloc((size_t)adj->n * sizeof(double));
    if (!dist) return false;
    bool result = !bellman_ford(adj, 0, dist);
    free(dist);
    return result;
}

apsp_path_t *apsp_extract_path(const apsp_graph_t *g, int32_t u, int32_t v) {
    if (!g || u < 0 || u >= g->n || v < 0 || v >= g->n) return NULL;
    apsp_path_t *p = malloc(sizeof(apsp_path_t));
    if (!p) return NULL;
    p->vertices = malloc((size_t)g->n * sizeof(int32_t));
    if (!p->vertices) { free(p); return NULL; }
    p->length = 0; p->total_weight = 0.0;
    if (g->dist[u * g->n + v] >= APSP_INF / 2.0) return p;
    int32_t cur = u;
    p->vertices[p->length++] = cur;
    while (cur != v && p->length < g->n) {
        int32_t nxt = g->next[cur * g->n + v];
        if (nxt < 0) break;
        p->vertices[p->length++] = nxt;
        cur = nxt;
    }
    p->total_weight = g->dist[u * g->n + v];
    return p;
}

double apsp_graph_diameter(const apsp_graph_t *g) {
    if (!g || g->n == 0) return 0.0;
    double max_d = 0.0;
    for (int32_t i = 0; i < g->n; i++)
        for (int32_t j = 0; j < g->n; j++)
            if (g->dist[i * g->n + j] < APSP_INF / 2.0 && g->dist[i * g->n + j] > max_d)
                max_d = g->dist[i * g->n + j];
    return max_d;
}

apsp_matrix_t *apsp_min_plus_multiply(const apsp_matrix_t *A, const apsp_matrix_t *B) {
    if (!A || !B || A->cols != B->rows) return NULL;
    int32_t r = A->rows, c = B->cols, inner = A->cols;
    apsp_matrix_t *C = malloc(sizeof(apsp_matrix_t));
    if (!C) return NULL;
    C->rows = r; C->cols = c;
    C->data = malloc((size_t)(r * c) * sizeof(double));
    if (!C->data) { free(C); return NULL; }
    for (int32_t i = 0; i < r; i++) {
        for (int32_t j = 0; j < c; j++) {
            double best = APSP_INF;
            for (int32_t k = 0; k < inner; k++) {
                double a = A->data[i * inner + k];
                double b = B->data[k * c + j];
                if (a < APSP_INF / 2.0 && b < APSP_INF / 2.0) {
                    double s = a + b;
                    if (s < best) best = s;
                }
            }
            C->data[i * c + j] = best;
        }
    }
    return C;
}

/* === Additional APSP Utilities and Analysis ========================== */

void apsp_print_routing_table(const apsp_graph_t *g, int32_t src) {
    if (!g || src < 0 || src >= g->n) return;
    printf("Routing table for vertex %d:\n", src);
    printf("  %-6s %-12s %-8s\n", "Dest", "Distance", "NextHop");
    for (int32_t v = 0; v < g->n; v++) {
        double d = g->dist[src * g->n + v];
        int32_t nh = g->next[src * g->n + v];
        if (d >= APSP_INF / 2.0) printf("  %-6d %-12s %-8s\n", v, "INF", "-");
        else printf("  %-6d %-12.4f %-8d\n", v, d, nh);
    }
}

double apsp_average_shortest_path(const apsp_graph_t *g) {
    if (!g || g->n < 2) return 0.0;
    double sum = 0.0; int32_t count = 0;
    for (int32_t i = 0; i < g->n; i++) {
        for (int32_t j = 0; j < g->n; j++) {
            if (i != j && g->dist[i * g->n + j] < APSP_INF / 2.0) {
                sum += g->dist[i * g->n + j]; count++;
            }
        }
    }
    return count > 0 ? sum / count : 0.0;
}

int32_t apsp_reachable_pairs(const apsp_graph_t *g) {
    if (!g) return 0;
    int32_t count = 0;
    for (int32_t i = 0; i < g->n; i++)
        for (int32_t j = 0; j < g->n; j++)
            if (i != j && g->dist[i * g->n + j] < APSP_INF / 2.0) count++;
    return count;
}

double apsp_eccentricity(const apsp_graph_t *g, int32_t v) {
    if (!g || v < 0 || v >= g->n) return 0.0;
    double ecc = 0.0;
    for (int32_t u = 0; u < g->n; u++)
        if (g->dist[v * g->n + u] < APSP_INF / 2.0 &&
            g->dist[v * g->n + u] > ecc)
            ecc = g->dist[v * g->n + u];
    return ecc;
}

double apsp_graph_radius(const apsp_graph_t *g) {
    if (!g || g->n == 0) return APSP_INF;
    double r = APSP_INF;
    for (int32_t v = 0; v < g->n; v++) {
        double ecc = apsp_eccentricity(g, v);
        if (ecc > 0.0 && ecc < r) r = ecc;
    }
    return r;
}

int32_t apsp_graph_center(const apsp_graph_t *g) {
    if (!g || g->n == 0) return -1;
    double min_ecc = APSP_INF; int32_t center = -1;
    for (int32_t v = 0; v < g->n; v++) {
        double ecc = apsp_eccentricity(g, v);
        if (ecc > 0.0 && ecc < min_ecc) { min_ecc = ecc; center = v; }
    }
    return center;
}

bool apsp_is_connected(const apsp_graph_t *g) {
    if (!g) return false;
    for (int32_t i = 0; i < g->n; i++)
        for (int32_t j = 0; j < g->n; j++)
            if (i != j && g->dist[i * g->n + j] >= APSP_INF / 2.0)
                return false;
    return true;
}

double apsp_network_routing_cost(const apsp_graph_t *g) {
    return apsp_average_shortest_path(g) * g->n;
}

double apsp_traffic_bottleneck(const apsp_graph_t *g) {
    return apsp_graph_diameter(g);
}

void apsp_print_adjacency(const apsp_adjacency_t *adj) {
    if (!adj) { printf("NULL adjacency\n"); return; }
    printf("Adjacency matrix (%dx%d):\n", adj->n, adj->n);
    int32_t display = adj->n < 10 ? adj->n : 10;
    for (int32_t i = 0; i < display; i++) {
        for (int32_t j = 0; j < display; j++) {
            double w = adj->weights[i * adj->n + j];
            if (w >= APSP_INF / 2.0) printf("  INF");
            else printf(" %4.1f", w);
        }
        printf("\n");
    }
    if (adj->n > 10) printf("  ... (truncated)\n");
}

int32_t apsp_edge_count(const apsp_adjacency_t *adj) {
    if (!adj) return 0;
    int32_t count = 0;
    for (int32_t i = 0; i < adj->n; i++)
        for (int32_t j = 0; j < adj->n; j++)
            if (i != j && adj->weights[i * adj->n + j] < APSP_INF / 2.0)
                count++;
    return count;
}

void apsp_print_statistics(const apsp_graph_t *g, const apsp_adjacency_t *adj) {
    printf("=== Graph Statistics ===\n");
    printf("Vertices: %d\n", g ? g->n : (adj ? adj->n : 0));
    printf("Edges:    %d\n", adj ? apsp_edge_count(adj) : 0);
    if (g) {
        printf("Diameter: %.4f\n", apsp_graph_diameter(g));
        printf("Radius:   %.4f\n", apsp_graph_radius(g));
        int32_t c = apsp_graph_center(g);
        printf("Center:   %s\n", c >= 0 ? "found" : "none");
        printf("Avg path: %.4f\n", apsp_average_shortest_path(g));
        printf("Connected: %s\n", apsp_is_connected(g) ? "yes" : "no");
        printf("Neg cycle: %s\n", g->has_negative_cycle ? "yes" : "no");
    }
}

bool apsp_verify_distances(const apsp_graph_t *g, const apsp_adjacency_t *adj) {
    if (!g || !adj || g->n != adj->n) return false;
    int32_t n = g->n;
    for (int32_t u = 0; u < n; u++) {
        for (int32_t v = 0; v < n; v++) {
            double direct = adj->weights[u * n + v];
            if (direct < APSP_INF / 2.0 && g->dist[u * n + v] > direct + 1e-12)
                return false;
            if (g->dist[u * n + v] < 0.0 && !g->has_negative_cycle)
                return false;
        }
    }
    return true;
}

apsp_graph_t *apsp_bellman_ford_all(const apsp_adjacency_t *adj) {
    if (!adj) return NULL;
    int32_t n = adj->n;
    apsp_graph_t *g = apsp_graph_create(n);
    if (!g) return NULL;
    double *dist = malloc((size_t)n * sizeof(double));
    if (!dist) { apsp_graph_destroy(g); return NULL; }
    for (int32_t s = 0; s < n; s++) {
        bool ok = bellman_ford(adj, s, dist);
        if (!ok) g->has_negative_cycle = true;
        for (int32_t t = 0; t < n; t++) g->dist[s * n + t] = dist[t];
    }
    free(dist); return g;
}
void apsp_print_path(const apsp_graph_t *g, int32_t u, int32_t v) {
    apsp_path_t *p = apsp_extract_path(g, u, v);
    if (!p) { printf("No path\n"); return; }
    if (p->length == 0) { printf("%d->%d: unreachable\n", u, v); free(p->vertices); free(p); return; }
    printf("%d->%d (dist=%.4f): ", u, v, p->total_weight);
    for (int32_t i = 0; i < p->length; i++) printf("%s%d", i>0?"->":"", p->vertices[i]);
    printf("\n"); free(p->vertices); free(p);
}
double apsp_betweenness_centrality(const apsp_graph_t *g, int32_t v) {
    if (!g || v<0 || v>=g->n) return 0.0;
    double bc = 0.0;
    for (int32_t s=0; s<g->n; s++) {
        for (int32_t t=0; t<g->n; t++) {
            if (s==t||s==v||t==v) continue;
            if (g->dist[s*g->n+v] >= APSP_INF/2.0) continue;
            if (g->dist[v*g->n+t] >= APSP_INF/2.0) continue;
            if (fabs(g->dist[s*g->n+v]+g->dist[v*g->n+t]-g->dist[s*g->n+t])<1e-12) bc+=1.0;
        }
    }
    return bc;
}
int32_t apsp_vertex_with_highest_betweenness(const apsp_graph_t *g) {
    if (!g||g->n==0) return -1;
    int32_t best=0; double best_bc=0.0;
    for (int32_t v=0; v<g->n; v++) {
        double bc=apsp_betweenness_centrality(g,v);
        if (bc>best_bc) { best_bc=bc; best=v; }
    }
    return best;
}
bool apsp_is_metric(const apsp_graph_t *g) {
    if (!g) return false;
    for (int32_t i=0; i<g->n; i++)
        for (int32_t j=0; j<g->n; j++)
            for (int32_t k=0; k<g->n; k++) {
                double dij=g->dist[i*g->n+j], dik=g->dist[i*g->n+k], dkj=g->dist[k*g->n+j];
                if (dij<APSP_INF/2.0&&dik<APSP_INF/2.0&&dkj<APSP_INF/2.0&&dij>dik+dkj+1e-12) return false;
            }
    return true;
}

double apsp_closeness_centrality(const apsp_graph_t *g, int32_t v) {
    if (!g || v<0 || v>=g->n) return 0.0;
    double sum=0.0; int32_t count=0;
    for (int32_t u=0; u<g->n; u++) {
        if (u!=v && g->dist[v*g->n+u] < APSP_INF/2.0) { sum+=g->dist[v*g->n+u]; count++; }
    }
    return count>0 ? (double)(count)/(sum*g->n) : 0.0;
}
int32_t apsp_count_triangles(const apsp_adjacency_t *adj) {
    if (!adj) return 0;
    int32_t n=adj->n, count=0;
    for (int32_t i=0; i<n; i++)
        for (int32_t j=i+1; j<n; j++)
            for (int32_t k=j+1; k<n; k++)
                if (adj->weights[i*n+j]<APSP_INF/2.0 && adj->weights[j*n+k]<APSP_INF/2.0 && adj->weights[k*n+i]<APSP_INF/2.0)
                    count++;
    return count;
}
bool apsp_has_negative_triangle(const apsp_adjacency_t *adj) {
    if (!adj) return false;
    int32_t n=adj->n;
    for (int32_t i=0; i<n; i++)
        for (int32_t j=0; j<n; j++)
            for (int32_t k=0; k<n; k++)
                if (adj->weights[i*n+j]<APSP_INF/2.0 && adj->weights[j*n+k]<APSP_INF/2.0 && adj->weights[k*n+i]<APSP_INF/2.0)
                    if (adj->weights[i*n+j]+adj->weights[j*n+k]+adj->weights[k*n+i] < -1e-12)
                        return true;
    return false;
}

void apsp_self_test(void) {
    printf("=== APSP Self-Test ===\n");
    int32_t n = 5;
    apsp_adjacency_t *adj = apsp_adjacency_create(n);
    apsp_adjacency_set_edge(adj, 0, 1, 3.0);
    apsp_adjacency_set_edge(adj, 0, 2, 8.0);
    apsp_adjacency_set_edge(adj, 1, 2, 4.0);
    apsp_adjacency_set_edge(adj, 1, 3, 1.0);
    apsp_adjacency_set_edge(adj, 2, 4, 2.0);
    apsp_adjacency_set_edge(adj, 3, 4, 7.0);
    apsp_adjacency_set_edge(adj, 1, 4, 6.0);
    apsp_graph_t *fw = apsp_floyd_warshall(adj);
    printf("Floyd-Warshall: dist(0,4)=%.1f (expected 9.0) %s\n",
           fw->dist[0*n+4], fabs(fw->dist[0*n+4]-9.0)<1e-6?"PASS":"FAIL");
    apsp_graph_t *dij = apsp_dijkstra_all(adj);
    printf("Dijkstra-All:   dist(0,4)=%.1f %s\n",
           dij->dist[0*n+4], fabs(dij->dist[0*n+4]-9.0)<1e-6?"PASS":"FAIL");
    printf("Diameter: %.1f\n", apsp_graph_diameter(fw));
    printf("Connected: %s\n", apsp_is_connected(fw)?"yes":"no");
    apsp_graph_destroy(fw); apsp_graph_destroy(dij);
    apsp_adjacency_destroy(adj);
}
double apsp_harmonic_centrality(const apsp_graph_t *g, int32_t v) {
    if (!g || v<0 || v>=g->n) return 0.0;
    double sum=0.0;
    for (int32_t u=0; u<g->n; u++)
        if (u!=v && g->dist[v*g->n+u] < APSP_INF/2.0 && g->dist[v*g->n+u] > 1e-12)
            sum += 1.0 / g->dist[v*g->n+u];
    return sum;
}

void apsp_explain_conjecture(void) {
    apsp_conjecture_t c = apsp_conjecture_status();
    printf("=== APSP Conjecture ===\n");
    printf("Statement: %s\n", c.statement);
    printf("History: %s\n", c.history);
}
void apsp_compare_algorithms(const apsp_adjacency_t *adj) {
    if (!adj) return;
    printf("=== APSP Algorithm Comparison (n=%d) ===\n", adj->n);
    printf("Algorithm          | Time           | Space      | Notes\n");
    printf("-------------------+----------------+------------+------\n");
    printf("Floyd-Warshall     | O(n^3)         | O(n^2)     | Dense graphs\n");
    printf("Johnson            | O(nm + n^2 log n)| O(n^2)  | Sparse graphs\n");
    printf("Dijkstra-All       | O(n(m + n log n))| O(n^2)  | No neg edges\n");
    printf("Bellman-Ford-All   | O(n^2 m)       | O(n^2)     | Neg edges OK\n");
    printf("Repeated Squaring  | O(M(n) log n)  | O(n^2)     | Min-plus mult\n");
    int32_t m = apsp_edge_count(adj);
    printf("\nThis graph: n=%d, m=%d edges\n", adj->n, m);
    printf("Floyd-Warshall ops: %d\n", adj->n * adj->n * adj->n);
    printf("Johnson ops:        %d (BF) + %d (Dijkstra)\n",
           adj->n * m, adj->n * (m + adj->n * adj->n));
}

/* End of apsp.c - All APSP algorithms implemented */
