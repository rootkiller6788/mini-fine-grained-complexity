/**
 * kclique_types.h — Core type definitions for k-Clique hardness module
 *
 * Defines graph representations, k-clique problem variants, and ETH/SETH
 * configuration types used throughout the fine-grained complexity analysis.
 *
 * Reference: Downey & Fellows, "Parameterized Complexity" (1999), Ch. 2, 7
 * Reference: Chen et al., "Tight Lower Bounds for Parameterized k-Clique" (2006)
 * Reference: Impagliazzo & Paturi, "Complexity of k-SAT" (1999/2001)
 *
 * Knowledge: L1 Definitions — k-Clique decision/search/counting problem types
 * Knowledge: L3 Mathematical Structures — Graph adjacency matrix, edge sets
 */

#ifndef KCLIQUE_TYPES_H
#define KCLIQUE_TYPES_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ================================================================
 * Graph representation (adjacency matrix)
 *
 * For n-vertex graph, adjacency is stored as a compact
 * matrix representation: adj[i * n + j] = 1 if edge (i,j) exists.
 * This supports O(1) adjacency queries essential for k-Clique.
 * ================================================================ */

typedef struct {
    uint8_t  *adj;        /* adjacency matrix, row-major: adj[i*n + j] */
    int32_t   n;          /* number of vertices */
    int32_t   m;          /* number of edges (cached for O(1) access) */
    bool      directed;   /* true if directed graph */
    bool      simple;     /* true if no self-loops or multi-edges */
} graph_t;

/**
 * k-Clique witness: a specific set of k vertices forming a clique.
 * vertices[0..k-1] are vertex indices in ascending order.
 */
typedef struct {
    int32_t  *vertices;   /* array of k vertex indices forming the clique */
    int32_t   k;          /* clique size */
    bool      found;      /* true if a clique was found */
} clique_t;

/**
 * Enumeration of clique variants studied in parameterized complexity.
 */
typedef enum {
    CLIQUETYPE_DECISION,   /* Does G contain a k-clique? (decision) */
    CLIQUETYPE_SEARCH,     /* Find one k-clique if it exists (search) */
    CLIQUETYPE_COUNT,      /* Count all k-cliques in G (counting) */
    CLIQUETYPE_LIST,       /* Enumerate all k-cliques (enumeration) */
    CLIQUETYPE_MAXIMAL,    /* Find a maximal clique (inclusion-maximal) */
    CLIQUETYPE_MAXIMUM,    /* Find a maximum clique (largest cardinality) */
    CLIQUETYPE_PARAMETERIZED /* Parameterized version with FPT verification */
} cliquetype_t;

/**
 * Bounded-degree clique parameters.
 * Monien's result: k-Clique in bounded-degree graphs is FPT
 * with parameter k (theory: O(k^2 d^k n) algorithm).
 */
typedef struct {
    int32_t   max_degree;   /* Delta(G), maximum vertex degree */
    bool      is_bounded;   /* true if degree is bounded by a constant */
    int32_t   degeneracy;   /* graph degeneracy (d-core number) */
    int32_t   arboricity;   /* graph arboricity alpha(G) */
} bounded_degree_t;

/**
 * ETH (Exponential Time Hypothesis) parameters.
 * ETH states: 3SAT cannot be solved in subexponential time 2^{o(n)}.
 * Consequence: k-Clique requires n^{Omega(k)} time.
 */
typedef struct {
    double    s_k;          /* liminf (log T(3SAT,n))/n  as n->inf */
    double    lower_bound;  /* lower bound exponent for k-Clique */
    int32_t   clique_parameter; /* k parameter for the clique problem */
    bool      eth_assumed;  /* whether ETH is assumed true */
} eth_params_t;

/**
 * SETH (Strong Exponential Time Hypothesis) parameters.
 * SETH: For any eps>0, there exists k such that kSAT requires 2^{(1-eps)n} time.
 * Consequence: k-Clique requires n^{k-o(1)} time for any constant k.
 */
typedef struct {
    double    cn;           /* constant c such that T(kSAT) >= 2^{c n} */
    double    epsilon;      /* eps parameter: T(kSAT) >= 2^{(1-eps)n} */
    int32_t   sat_k;        /* kSAT parameter triggering SETH consequence */
    int32_t   clique_k;     /* k parameter for clique lower bound */
    double    implied_lower; /* implied n^{implied_lower} lower bound */
    bool      seth_assumed; /* whether SETH is assumed true */
} seth_params_t;

/**
 * Conditional lower bound result: encodes a theorem of the form
 * "If conjecture X is true, then problem P requires Omega(n^alpha) time."
 */
typedef struct {
    char     *conjecture;       /* name of the hardness conjecture */
    char     *problem;          /* target problem */
    double    time_exponent;    /* alpha: n^alpha lower bound, or INFINITY */
    bool      exponential;      /* true if lower bound is 2^{Omega(n)} */
    double    conditional_constant; /* leading constant factor */
    bool      tight;            /* whether the bound is known to be tight */
} conditional_lower_bound_t;

/**
 * W[1]-hardness certificate for k-Clique.
 * Downey & Fellows (1995): k-Clique is W[1]-complete.
 * This struct captures the FPT reduction parameters.
 */
typedef struct {
    int32_t   original_k;       /* parameter of the W[1]-hard source problem */
    int32_t   clique_k;         /* k value in the constructed k-Clique instance */
    int32_t   instance_size;    /* n: number of vertices in constructed graph */
    double    reduction_time;   /* f(k) * n^{O(1)}: FPT reduction runtime */
    bool      is_fpt_reduction; /* true if the reduction is an FPT reduction */
} w1_hardness_cert_t;

/**
 * Graph generation parameters for Erdos-Renyi model G(n,p).
 * Useful for testing and empirical analysis of clique algorithms.
 */
typedef struct {
    int32_t   n;              /* number of vertices */
    double    p;              /* edge probability */
    int64_t   seed;           /* random seed for reproducibility */
    bool      force_connected; /* ensure the graph is connected */
    bool      include_clique;  /* plant a hidden clique for testing */
    int32_t   planted_k;       /* size of planted clique */
} er_graph_params_t;

/**
 * Color-coding algorithmic state.
 * Used by Alon-Yuster-Zwick (1995) color-coding FPT algorithm.
 * Runtime: O(2^k * n^2) for k-Clique.
 */
typedef struct {
    int32_t  *colors;          /* colors[v] in {0,...,k-1} for each vertex */
    int32_t   n;               /* number of vertices */
    int32_t   k;               /* target clique size */
    int32_t   num_iterations;  /* number of random colorings attempted */
    int64_t   seed;            /* RNG seed */
    bool      use_explicit_hash; /* use explicit k-perfect hash family */
} color_coding_t;

/**
 * Branch-and-bound state for maximum clique search.
 * Tomita et al. (2003/2016): efficient branch-and-bound for MCP.
 */
typedef struct {
    int32_t  *current_solution;  /* vertices in current candidate */
    int32_t   current_size;      /* |current_solution| */
    int32_t  *best_solution;     /* best clique found so far */
    int32_t   best_size;         /* |best_solution| */
    int64_t   nodes_explored;    /* search tree nodes evaluated */
    int32_t   upper_bound_limit; /* pruning threshold */
    bool      use_approximate_coloring; /* use greedy coloring for bound */
} branch_and_bound_t;

/**
 * Ramsey numbers R(s,t): the smallest n such that any n-vertex graph
 * contains either a clique of size s or an independent set of size t.
 * Ramsey (1930), Erdos-Szekeres (1935).
 */
typedef struct {
    int32_t   s;              /* clique size parameter */
    int32_t   t;              /* independent set size parameter */
    int64_t   lower_bound;    /* known lower bound on R(s,t) */
    int64_t   upper_bound;    /* known upper bound on R(s,t) */
    bool      exact_known;    /* true if R(s,t) is exactly known */
} ramsey_number_t;

/**
 * Circuit representation for W[1]-hardness analysis.
 * W[1] is defined via Boolean circuits of bounded weft.
 * A circuit has weft t if the maximum number of gates with
 * fan-in > 2 on any input-to-output path is t.
 */
typedef struct {
    int32_t   num_gates;      /* total gates in circuit */
    int32_t   num_inputs;     /* number of input variables */
    int32_t   weft;           /* weft (depth of large fan-in gates) */
    int32_t   depth;          /* total circuit depth */
    bool      is_monotone;    /* true if all gates are AND/OR (no NOT) */
    bool      is_bounded_fanin; /* true if all gates have bounded fan-in */
} circuit_weft_t;

/**
 * The k-Clique problem configuration bundles problem type,
 * graph, and algorithm parameters for modular function dispatch.
 */
typedef struct {
    graph_t          graph;          /* input graph */
    int32_t          k;              /* clique size parameter */
    cliquetype_t     problem_type;   /* decision/search/count */
    eth_params_t     eth;            /* ETH parameters */
    seth_params_t    seth;           /* SETH parameters */
    bool             use_color_coding; /* enable color-coding algorithm */
    bool             use_branch_bound; /* enable branch-and-bound */
    int32_t          report_interval;  /* progress report frequency */
} kclique_config_t;

/* ================================================================
 * Parameterized complexity: Fixed-Parameter Tractability
 *
 * A parameterized problem L subseteq Sigma* x N is FPT if there exists
 * an algorithm deciding (x,k) in L in time f(k)*|x|^{O(1)} for some
 * computable f.
 * ================================================================ */

typedef struct {
    int32_t   kernel_size;         /* size of kernel (reduced instance) */
    double    kernelization_time;  /* time to compute kernel */
    int32_t   original_size;       /* |x|: original instance size */
    int32_t   parameter;           /* k: parameter value */
    bool      polynomial_kernel;   /* true if kernel size is poly(k) */
    bool      admits_kernelization; /* whether problem admits a kernel */
} fpt_kernel_t;

/**
 * FPT reduction from problem A to problem B.
 * An FPT reduction (also called parameterized reduction) maps
 * (x,k) -> (x',k') such that:
 * 1. (x,k) in A iff (x',k') in B
 * 2. k' <= g(k) for some computable g
 * 3. The reduction runs in FPT time f(k)*|x|^{O(1)}
 */
typedef struct {
    char     *source_problem;     /* problem A */
    char     *target_problem;     /* problem B */
    double    runtime_exponent;   /* |x|^c factor exponent */
    int32_t   param_blowup;       /* g(k) bound on parameter growth */
    bool      is_fpt;             /* true if reduction is FPT */
} fpt_reduction_t;

/* ================================================================
 * Configuration validation flags
 * ================================================================ */

#define GRAPH_MAX_VERTICES   (1 << 20)   /* max n = 1,048,576 */
#define GRAPH_MAX_K          128         /* max k for exact algorithms */
#define ETH_DEFAULT_SK       1.0         /* s_k >= 1 (ETH: no 2^{o(n)} for 3SAT) */
#define SETH_DEFAULT_CN      1.0         /* c_n -> 1 as n -> inf (SETH) */

#ifdef __cplusplus
}
#endif



/* ================================================================
 * Utility functions shared across modules
 * ================================================================ */

/** Check if a kernel size qualifies as polynomial (size <= k^10). */
bool is_polynomial_kernel(int32_t kernel_size, int32_t k);

/** Compute lower bound on kernel size for k-Clique under ETH. */
int64_t kclique_kernel_lower_bound(int32_t k);

/** Analyze parameterized complexity for k-Clique on a graph class. */
typedef enum {
    PARAM_RESULT_FPT,
    PARAM_RESULT_W1_HARD,
    PARAM_RESULT_XP,
    PARAM_RESULT_PARA_NP,
    PARAM_RESULT_UNKNOWN
} param_result_t;

param_result_t kclique_param_complexity(int32_t max_degree,
                                         bool is_planar,
                                         bool is_bipartite,
                                         bool has_bounded_treewidth);

#endif /* KCLIQUE_TYPES_H */
