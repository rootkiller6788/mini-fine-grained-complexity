/* Example: 3SUM Equivalence -- 3SUM, Collinearity, Applications
 *
 * Demonstrates the 3SUM equivalence class:
 *   1. Standard 3SUM (sort + two-pointer)
 *   2. Collinearity detection (3SUM-hard)
 *   3. Mars rover obstacle avoidance (NASA application)
 *   4. Right triangle detection
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "threesum.h"
#include <string.h>

int main(void) {
    printf("=== Example: 3SUM Equivalence Class ===\n\n");

    /* --- 3SUM --- */
    printf("--- 3SUM Problem ---\n");
    threesum_instance_t *ts = threesum_create(8);
    int64_t vals[] = {-25, -10, -7, -3, 2, 4, 8, 10};
    for (int32_t i = 0; i < 8; i++) ts->values[i] = vals[i];

    printf("Set: { -25, -10, -7, -3, 2, 4, 8, 10 }\n");
    bool found = threesum_quadratic(ts);
    printf("3SUM: %s\n", found ? "YES" : "NO");
    if (found) printf("  Solution: %lld + %lld + %lld = 0\n",
        (long long)ts->values[ts->idx_a],
        (long long)ts->values[ts->idx_b],
        (long long)ts->values[ts->idx_c]);
    threesum_free(ts);

    /* --- 3SUM with no solution --- */
    printf("\nSet: { 1, 2, 3, 4, 5 }\n");
    threesum_instance_t *ts2 = threesum_create(5);
    ts2->values[0] = 1; ts2->values[1] = 2; ts2->values[2] = 3;
    ts2->values[3] = 4; ts2->values[4] = 5;
    found = threesum_quadratic(ts2);
    printf("3SUM: %s\n", found ? "YES (unexpected!)" : "NO (as expected)");
    threesum_free(ts2);

    /* --- Collinearity Detection --- */
    printf("\n--- Collinearity Detection (3SUM-hard) ---\n");
    collinearity_instance_t *ci = collinearity_create(5);
    double cx[] = {0.0, 1.0, 2.0, 3.0, 5.0};
    double cy[] = {0.0, 1.0, 2.0, 3.0, 5.0};
    memcpy(ci->x, cx, 5 * sizeof(double));
    memcpy(ci->y, cy, 5 * sizeof(double));
    bool col = collinearity_via_3sum(ci);
    printf("Collinear triple: %s\n", col ? "YES" : "NO");
    if (col) printf("  Points: (%g,%g), (%g,%g), (%g,%g)\n",
        cx[ci->i1], cy[ci->i1], cx[ci->i2], cy[ci->i2], cx[ci->i3], cy[ci->i3]);
    collinearity_free(ci);

    /* --- Mars Rover Obstacle Check --- */
    printf("\n--- Mars Rover Obstacle Planning (NASA) ---\n");
    double ox[] = {0.0, 2.0, 4.0, 6.0, 1.0};
    double oy[] = {0.0, 2.0, 4.0, 6.0, 0.5};
    bool blocked = mars_rover_collinearity_obstacle_check(ox, oy, 5);
    printf("Collinear obstacles blocking path: %s\n", blocked ? "YES" : "NO");

    /* --- Right Triangle Detection --- */
    printf("\n--- Right Triangle Detection ---\n");
    double rx[] = {0.0, 3.0, 0.0, 1.0};
    double ry[] = {0.0, 0.0, 4.0, 1.0};
    int32_t i1, i2, i3;
    bool rt = right_triangle_detection(rx, ry, 4, &i1, &i2, &i3);
    printf("Right triangle: %s\n", rt ? "YES" : "NO");
    if (rt) printf("  Vertices: (%g,%g), (%g,%g), (%g,%g)\n",
        rx[i1], ry[i1], rx[i2], ry[i2], rx[i3], ry[i3]);

    /* --- Status Report --- */
    printf("\n--- 3SUM Status ---\n");
    threesum_status_report();

    printf("\n=== Example Complete ===\n");
    return 0;
}