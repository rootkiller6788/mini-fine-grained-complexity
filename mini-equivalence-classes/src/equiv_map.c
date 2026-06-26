/* ============================================================================
 * equiv_map.c -- Equivalence Class Mapping and Problem Registration
 *
 * Initializes the known fine-grained equivalence classes, registers
 * canonical problems, and establishes cross-class reduction relationships.
 *
 * This file serves as the "database" of known fine-grained results:
 *   - Subcubic class members and their relationships
 *   - Subquadratic class members and their relationships
 *   - 3SUM class members and their relationships
 *   - Cross-class implications (e.g., SETH => OV => Edit Distance)
 *
 * Knowledge Coverage:
 *   L1: Problem registration for all canonical problems
 *   L2: Fine-grained reduction registration
 *   L4: Key theorems as reduction edges
 *   L6: All canonical problems initialized
 *   L7: Practical applications mapped to theoretical problems
 *   L8: Cross-class landscape analysis
 *   L9: Current open problems and conjectures documented
 *
 * References:
 *   Williams & Williams (2013), STOC -- Subcubic Equivalences
 *   Williams (2005), ICALP -- SETH => OV
 *   Backurs & Indyk (2016), STOC -- OV => Edit Distance
 *   Bringmann (2014), FOCS -- OV => Frechet
 *   Abboud, Backurs, Williams (2015), FOCS -- OV => LCS, DTW
 * ============================================================================ */

#include "equiv_classes.h"
#include "subcubic.h"
#include "subquadratic.h"
#include "threesum.h"
#include "fine_grained_reduction.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ---- Problem ID registry ---- */

static problem_id_t PID_APSP;
static problem_id_t PID_NEGTRI;
static problem_id_t PID_MINPLUS;
static problem_id_t PID_RADIUS;
static problem_id_t PID_MEDIAN;
static problem_id_t PID_DIAMETER;
static problem_id_t PID_OV;
static problem_id_t PID_EDITDIST;
static problem_id_t PID_LCS;
static problem_id_t PID_FRECHET;
static problem_id_t PID_DTW;
static problem_id_t PID_3SUM;
static problem_id_t PID_COLLINEAR;

/* ---- Initialize all known equivalence classes ---- */

static void init_equivalence_classes(void) {
    /* Subcubic class: canonical = APSP, threshold = 3.0 */
    equiv_class_create(EQUIV_CLASS_SUBCUBIC, "APSP", 3.0, 3.0);
    equiv_class_create(EQUIV_CLASS_APSP, "APSP", 3.0, 3.0);
    equiv_class_create(EQUIV_CLASS_BMM, "Boolean Matrix Multiplication", 3.0, 3.0);
    equiv_class_create(EQUIV_CLASS_DIAMETER, "Graph Diameter", 3.0, 3.0);

    /* Subquadratic class: canonical = OV, threshold = 2.0 */
    equiv_class_create(EQUIV_CLASS_SUBQUADRATIC, "Orthogonal Vectors", 2.0, 2.0);

    /* 3SUM class: canonical = 3SUM, threshold = 2.0 */
    equiv_class_create(EQUIV_CLASS_3SUM, "3SUM", 2.0, 2.0);
    equiv_class_create(EQUIV_CLASS_COLLINEARITY, "Collinearity", 2.0, 2.0);

    /* Exponential class: canonical = CNF-SAT */
    equiv_class_create(EQUIV_CLASS_CNFSAT, "k-SAT", 1.0, 1.0);

    /* All-Pairs LCS */
    equiv_class_create(EQUIV_CLASS_ALLPAIRS_LCS, "All-Pairs LCS", 2.5, 2.5);

    /* Hitting Set */
    equiv_class_create(EQUIV_CLASS_HITTING_SET, "Hitting Set", 2.0, 2.0);
}

/* ---- Register all canonical problems ---- */

static void register_all_problems(void) {
    PID_APSP     = register_problem("APSP", EQUIV_CLASS_SUBCUBIC, 3.0, 10);
    PID_NEGTRI   = register_problem("Negative Triangle", EQUIV_CLASS_SUBCUBIC, 3.0, 10);
    PID_MINPLUS  = register_problem("Min-Plus Product", EQUIV_CLASS_SUBCUBIC, 3.0, 10);
    PID_RADIUS   = register_problem("Graph Radius", EQUIV_CLASS_DIAMETER, 3.0, 10);
    PID_MEDIAN   = register_problem("Graph Median", EQUIV_CLASS_DIAMETER, 3.0, 10);
    PID_DIAMETER = register_problem("Graph Diameter", EQUIV_CLASS_DIAMETER, 3.0, 10);

    PID_OV       = register_problem("Orthogonal Vectors", EQUIV_CLASS_SUBQUADRATIC, 2.0, 10);
    PID_EDITDIST = register_problem("Edit Distance", EQUIV_CLASS_SUBQUADRATIC, 2.0, 10);
    PID_LCS      = register_problem("LCS", EQUIV_CLASS_SUBQUADRATIC, 2.0, 10);
    PID_FRECHET  = register_problem("Frechet Distance", EQUIV_CLASS_SUBQUADRATIC, 2.0, 10);
    PID_DTW      = register_problem("DTW", EQUIV_CLASS_SUBQUADRATIC, 2.0, 10);

    PID_3SUM     = register_problem("3SUM", EQUIV_CLASS_3SUM, 2.0, 10);
    PID_COLLINEAR = register_problem("Collinearity", EQUIV_CLASS_COLLINEARITY, 2.0, 10);
}

/* ---- Register all known reductions ---- */

static void register_all_reductions(void) {
    /* ===== Subcubic class: APSP <-> NegTriangle <-> Min-Plus ===== */

    /* APSP -> NegTriangle (trivial via distance check) */
    register_reduction(PID_APSP, PID_NEGTRI, FGR_SUBCUBIC, 0.0, 1,
        "Williams & Williams (2013): APSP via NegTriangle oracle");
    /* NegTriangle -> APSP (via min-plus product) */
    register_reduction(PID_NEGTRI, PID_APSP, FGR_SUBCUBIC, 3.0, 1,
        "Williams & Williams (2013): NegTriangle reduces to APSP");
    /* APSP -> Min-Plus (repeated squaring) */
    register_reduction(PID_APSP, PID_MINPLUS, FGR_SUBCUBIC, 0.0, 1,
        "Folklore: APSP via log n min-plus products");
    /* Min-Plus -> APSP (direct) */
    register_reduction(PID_MINPLUS, PID_APSP, FGR_SUBCUBIC, 0.0, 1,
        "Folklore: Min-Plus is a single APSP step");

    /* APSP -> Radius */
    register_reduction(PID_APSP, PID_RADIUS, FGR_SUBCUBIC, 2.0, 1,
        "Williams & Williams (2013): Radius <=_subcubic APSP");
    /* APSP -> Median */
    register_reduction(PID_APSP, PID_MEDIAN, FGR_SUBCUBIC, 2.0, 1,
        "Williams & Williams (2013): Median <=_subcubic APSP");
    /* APSP -> Diameter */
    register_reduction(PID_APSP, PID_DIAMETER, FGR_SUBCUBIC, 2.0, 1,
        "Williams & Williams (2013): Diameter <=_subcubic APSP");

    /* ===== Subquadratic class: OV <-> EditDistance <-> LCS <-> DTW ===== */

    /* OV -> Edit Distance (Backurs & Indyk 2016) */
    register_reduction(PID_OV, PID_EDITDIST, FGR_SUBQUADRATIC, 2.0, 1,
        "Backurs & Indyk (2016): OV <=_subquadratic Edit Distance");
    /* Edit Distance -> OV */
    register_reduction(PID_EDITDIST, PID_OV, FGR_SUBQUADRATIC, 2.0, 1,
        "Subquadratic equivalence: Edit Distance => OV");

    /* OV -> LCS (Abboud, Backurs, Williams 2015) */
    register_reduction(PID_OV, PID_LCS, FGR_SUBQUADRATIC, 2.0, 1,
        "Abboud, Backurs, Williams (2015): OV <=_subquadratic LCS");
    register_reduction(PID_LCS, PID_OV, FGR_SUBQUADRATIC, 2.0, 1,
        "Subquadratic equivalence: LCS => OV");

    /* OV -> Frechet (Bringmann 2014) */
    register_reduction(PID_OV, PID_FRECHET, FGR_SUBQUADRATIC, 2.0, 1,
        "Bringmann (2014): OV <=_subquadratic Frechet");
    register_reduction(PID_FRECHET, PID_OV, FGR_SUBQUADRATIC, 2.0, 1,
        "Subquadratic equivalence: Frechet => OV");

    /* OV -> DTW (Abboud, Backurs, Williams 2015) */
    register_reduction(PID_OV, PID_DTW, FGR_SUBQUADRATIC, 2.0, 1,
        "Abboud et al. (2015): OV <=_subquadratic DTW");
    register_reduction(PID_DTW, PID_OV, FGR_SUBQUADRATIC, 2.0, 1,
        "Subquadratic equivalence: DTW => OV");

    /* ===== 3SUM class ===== */
    register_reduction(PID_3SUM, PID_COLLINEAR, FGR_3SUM, 2.0, 1,
        "Gajentaan & Overmars (1995): 3SUM <=_subquadratic Collinearity");
    register_reduction(PID_COLLINEAR, PID_3SUM, FGR_3SUM, 2.0, 1,
        "3SUM equivalence: Collinearity => 3SUM");
}

/* ---- Exported: Initialize the entire equivalence class framework ---- */

void equiv_classes_init(void) {
    static bool initialized = false;
    if (initialized) return;
    init_equivalence_classes();
    register_all_problems();
    register_all_reductions();
    initialized = true;
}

/* ---- Convenience: Print the full fine-grained complexity map ---- */

void equiv_classes_print_map(void) {
    equiv_classes_init();
    print_equivalence_class_summary();
    printf("\n");
    printf("=== Fine-Grained Reduction Map ===\n\n");
    printf("Subcubic Equivalences (APSP ~ NegTriangle ~ MinPlus ~ Radius ~ Median ~ Diameter):\n");
    printf("  All pairs subcubic-equivalent (Williams & Williams 2013, STOC)\n\n");
    printf("Subquadratic Equivalences (OV ~ EditDistance ~ LCS ~ Frechet ~ DTW):\n");
    printf("  All pairs subquadratic-equivalent under SETH\n");
    printf("  Implied by: SETH => OV conjecture (Williams 2005)\n\n");
    printf("3SUM Equivalences (3SUM ~ Collinearity ~ PolygonFitting):\n");
    printf("  All pairs 3SUM-equivalent (Gajentaan & Overmars 1995)\n\n");
    printf("Cross-class relationships:\n");
    printf("  SETH => OV conjecture => All subquadratic quadratic lower bounds\n");
    printf("  3SUM conjecture is independent of SETH (weaker assumption)\n");
    printf("  Subcubic conjecture is independent of SETH\n");
}

/* ---- L9: Research Frontiers ---- */

void equiv_classes_print_open_problems(void) {
    printf("=== Open Problems in Fine-Grained Complexity ===\n\n");
    printf("1. Subcubic Conjecture:\n");
    printf("   Is APSP in O(n^{3-eps}) time? Known: O(n^3/log^2 n).\n");
    printf("   Would fast matrix multiplication (omega=2) refute it?\n\n");
    printf("2. OV Conjecture:\n");
    printf("   Is OV in O(n^{2-eps}) for d=omega(log n)?\n");
    printf("   Equivalent to SETH. Very strongly believed.\n\n");
    printf("3. 3SUM Conjecture:\n");
    printf("   Is 3SUM in O(n^{2-eps})? Known: O(n^2/polylog n).\n");
    printf("   True in decision tree model (Erickson 1999).\n\n");
    printf("4. APSP vs OV:\n");
    printf("   Does APSP hardness imply OV hardness?\n");
    printf("   Currently no known reduction either way.\n\n");
    printf("5. Combinatorial BMM:\n");
    printf("   Is there a combinatorial (non-algebraic) algorithm for\n");
    printf("   Boolean Matrix Multiplication in O(n^{3-eps})?\n");
    printf("   Conjectured to require n^{3-o(1)} time.\n\n");
    printf("6. Approximability:\n");
    printf("   Can Edit Distance be approximated in near-linear time?\n");
    printf("   Known: No O(n^{2-eps}) exact algorithm (under SETH),\n");
    printf("   but can be approximated within factor 1+eps in O(n log n)?\n");
    printf("   Open: Tight lower bounds for approximation.\n\n");
    printf("References:\n");
    printf("  Williams (2018), \"Some Open Problems in Fine-Grained Complexity\"\n");
    printf("  Rubinstein & Williams (2019), \"Open Problems\" column, SIGACT News\n");
}