/**
 * ex4_w1_hardness.c ? W[1]-Hardness of k-Clique Demonstration
 *
 * Demonstrates the 3SAT -> k-Clique FPT reduction that proves
 * k-Clique is W[1]-hard. Creates a small 3SAT instance, constructs
 * the reduction graph, finds a clique, and recovers the SAT assignment.
 *
 * Knowledge: L4 ? W[1]-completeness proof via reduction
 * Knowledge: L6 ? 3SAT, k-Clique canonical problems
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../include/kclique_types.h"
#include "../include/kclique_core.h"
#include "../include/kclique_reduction.h"
#include "../include/kclique_algorithm.h"

int main(void) {
    printf("=== W[1]-Hardness: 3SAT to k-Clique Reduction ===\n\n");
    printf("This demonstrates the canonical FPT reduction that proves\n");
    printf("k-Clique is W[1]-hard (Downey & Fellows, 1995).\n\n");

    /* Create a small satisfiable 3SAT instance:
       (x1 OR x2 OR x3) AND (NOT x1 OR x2 OR NOT x3) AND (x1 OR NOT x2 OR x3)
       Variables: 3, Clauses: 3
       Satisfying assignment: x1=T, x2=T, x3=T */
    int32_t nvars = 3, nclauses = 3;
    sat3_instance_t sat;
    sat.num_variables = nvars;
    sat.num_clauses = nclauses;
    sat.clause_literals = (int32_t *)malloc(3 * nclauses * sizeof(int32_t));

    /* Clause 1: x1 OR x2 OR x3  (all positive => vars 0,1,2 sign=true) */
    sat.clause_literals[0] = 0;   /* x1 positive */
    sat.clause_literals[1] = 1;   /* x2 positive */
    sat.clause_literals[2] = 2;   /* x3 positive */

    /* Clause 2: NOT x1 OR x2 OR NOT x3
       Positive literal: value >= 0 and positive
       Negative literal: value < 0, encode as -(var+1) */
    sat.clause_literals[3] = -1;  /* NOT x1 (var 0) */
    sat.clause_literals[4] = 1;   /* x2 positive */
    sat.clause_literals[5] = -3;  /* NOT x3 (var 2) */

    /* Clause 3: x1 OR NOT x2 OR x3 */
    sat.clause_literals[6] = 0;   /* x1 positive */
    sat.clause_literals[7] = -2;  /* NOT x2 (var 1) */
    sat.clause_literals[8] = 2;   /* x3 positive */

    printf("3SAT Instance:\n");
    printf("  Variables: %d\n", nvars);
    printf("  Clauses: %d\n", nclauses);
    printf("  (x1 OR x2 OR x3) AND (NOT x1 OR x2 OR NOT x3) AND (x1 OR NOT x2 OR x3)\n\n");

    /* Perform the reduction */
    int32_t k;
    graph_t *reduction_graph = reduce_3sat_to_kclique_fpt(&sat, &k);
    if (!reduction_graph) {
        printf("Reduction failed.\n");
        free(sat.clause_literals);
        return 1;
    }

    printf("Reduction constructed:\n");
    printf("  Graph vertices: %d (7 * %d = 21)\n", reduction_graph->n, nclauses);
    printf("  Target clique size: k = %d\n\n", k);

    /* Verify reduction properties */
    bool valid = verify_3sat_to_kclique_reduction(&sat, reduction_graph, k);
    printf("Reduction verification: %s\n\n", valid ? "PASSED" : "FAILED");

    /* Find a k-clique in the constructed graph */
    clique_t found_clique = {0};
    bool clique_found = color_coding_solve(reduction_graph, k, 200, 123, &found_clique);

    printf("k-Clique search: %s\n", clique_found ? "FOUND" : "NOT FOUND");

    if (clique_found) {
        printf("  Clique vertices: {");
        for (int i = 0; i < found_clique.k; i++) {
            printf("%d%s", found_clique.vertices[i],
                   i < found_clique.k - 1 ? ", " : "");
        }
        printf("}\n\n");

        /* Recover SAT assignment */
        bool *assignment = (bool *)calloc(nvars, sizeof(bool));
        if (recover_sat_assignment(&sat, &found_clique, assignment)) {
            printf("Recovered SAT assignment:\n");
            for (int v = 0; v < nvars; v++) {
                printf("  x%d = %s\n", v + 1, assignment[v] ? "TRUE" : "FALSE");
            }
        } else {
            printf("Failed to recover valid SAT assignment.\n");
        }
        free(assignment);
        clique_free(&found_clique);
    }

    printf("\nKey Insight:\n");
    printf("  This reduction proves k-Clique is W[1]-hard because:\n");
    printf("  1. 3SAT (parameterized by #clauses) is W[1]-hard\n");
    printf("  2. The reduction is parameter-preserving (k = #clauses)\n");
    printf("  3. The reduction runs in polynomial time\n");
    printf("  4. Therefore: if k-Clique in FPT, then W[1] = FPT\n");

    graph_destroy(reduction_graph);
    free(sat.clause_literals);
    printf("\nDone.\n");
    return 0;
}
