/* ============================================================================
 * main.c -- Mini Equivalence Classes: Entry Point
 *
 * Demonstrates the fine-grained equivalence class framework.
 * ============================================================================ */

#include <stdio.h>
#include "equiv_classes.h"
#include "subcubic.h"
#include "subquadratic.h"
#include "threesum.h"
#include "fine_grained_reduction.h"

/* Forward declarations from equiv_map.c */
extern void equiv_classes_init(void);
extern void equiv_classes_print_map(void);
extern void equiv_classes_print_open_problems(void);

int main(void) {
    printf("=== mini-equivalence-classes ===\n");
    printf("Fine-Grained Equivalence Classes of Problems\n\n");

    /* Initialize the framework */
    equiv_classes_init();

    /* Print the landscape */
    equiv_classes_print_map();
    printf("\n");

    /* Print status reports for each class */
    subcubic_status_report();
    printf("\n");
    subquadratic_status_report();
    printf("\n");
    threesum_status_report();
    printf("\n");

    /* Print open research problems */
    equiv_classes_print_open_problems();

    printf("\n=== Module Status: COMPLETE ===\n");
    return 0;
}