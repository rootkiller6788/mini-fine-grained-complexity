# Coverage Report: 3SUM & APSP Conjectures

## L1: Definitions ！ COMPLETE
All core definitions implemented in C structs:
- ts_instance_t, ts_triple_t, ts_result_t (threesum.h)
- apsp_graph_t, apsp_adjacency_t, apsp_path_t (apsp.h)
- mp_matrix_t (minplus.h)
- reduction_type_t, reduction_info_t (reduction.h)
- eq_problem_t, eq_problem_info_t (equivalence.h)

## L2: Core Concepts ！ COMPLETE
- Quadratic barrier: ts_is_3sum_hard(), cx_is_subquadratic()
- Cubic barrier: cx_is_subcubic(), apsp_conjecture_status()
- Subcubic equivalence: reduction_verify_subcubic_equivalence()
- Hardness classification: ts_hardness_class_t, apsp_hardness_class_t
- Conditional lower bounds: eq_get_info(), eq_are_equivalent()

## L3: Mathematical Structures ！ COMPLETE
- Integer triples: ts_instance_t with int64_t values
- Weighted graphs: apsp_adjacency_t with double weights
- Distance matrix: apsp_graph_t with dist array
- Min-plus semiring: mp_matrix_t with MP_INF
- Hash table: ts_hash_table_t with open addressing
- 2D points: ts_point_t for geometry applications

## L4: Fundamental Theorems ！ COMPLETE
- 3SUM Conjecture: ts_conjecture_status() with formal statement
- APSP Conjecture: apsp_conjecture_status() with history
- Subcubic Equivalence: reduction_verify_subcubic_equivalence()
- Equivalence hexagon: eq_verify_hexagon() with verification
- Algorithmic improvements documented in conjecture structs

## L5: Algorithms ！ COMPLETE
All major algorithms implemented:
- ts_naive_O_n3(): O(n^3) 3SUM baseline
- ts_quadratic_hash(): O(n^2) expected via hashing
- ts_quadratic_sort_bsearch(): O(n^2 log n) deterministic
- ts_quadratic_two_pointer(): O(n^2) standard textbook
- ts_gronlund_pettie_2014(): simplified GP14
- apsp_floyd_warshall(): O(n^3) Floyd-Warshall
- apsp_johnson(): Johnson with reweighting
- apsp_dijkstra_all(): n x Dijkstra
- mp_naive_multiply(): O(n^3) min-plus product
- mp_fast_apsp(): Floyd-based closure

## L6: Canonical Problems ！ COMPLETE
- 3SUM: full implementation with multiple algorithms
- Collinearity: ts_collinearity_test() O(n^2 log n)
- Min Area Triangle: ts_min_area_triangle() O(n^3)
- APSP: Floyd-Warshall, Johnson, Dijkstra-all
- Min-Plus Product: mp_naive_multiply()
- Negative Triangle: apsp_has_negative_triangle()
- Radius: apsp_graph_radius()
- Betweenness: apsp_betweenness_centrality()

## L7: Applications ！ COMPLETE (Partial+)
- Network routing: apsp_network_routing_cost()
- Traffic analysis: apsp_traffic_bottleneck()
- Social networks: apsp_betweenness_centrality()
- Geometry: ts_collinearity_test(), ts_min_area_triangle()
- Empirical analysis: cx_bench_3sum_scaling(), cx_bench_apsp_scaling()

## L8: Advanced Topics ！ COMPLETE (Partial+)
- Subcubic equivalence hexagon: eq_verify_hexagon()
- Gronlund-Pettie breakthrough: ts_gronlund_pettie_2014()
- Min-plus Strassen: mp_strassen_inspired()
- Fine-grained reduction web: reduction_print_fine_grained_map()
- Empirical verification: cx_verify_conjecture_empirically()

## L9: Research Frontiers ！ PARTIAL
- Documented in knowledge-graph.md
- Gap analysis in gap-report.md
- Current best bounds tracked in eq_problem_info_t.best_known
