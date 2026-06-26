# Knowledge Graph: mini-equivalence-classes

## L1: Definitions
| # | Definition | C Implementation | Status |
|---|-----------|-----------------|--------|
| 1 | Equivalence Class (`equiv_class_t`) | `include/equiv_classes.h` | Complete |
| 2 | Fine-Grained Reduction (`fg_reduction_type_t`) | `include/equiv_classes.h` | Complete |
| 3 | Problem Descriptor (`problem_descriptor_t`) | `include/equiv_classes.h` | Complete |
| 4 | Reduction Descriptor (`reduction_descriptor_t`) | `include/equiv_classes.h` | Complete |
| 5 | APSP Instance (`apsp_instance_t`) | `include/subcubic.h` | Complete |
| 6 | Negative Triangle (`neg_triangle_instance_t`) | `include/subcubic.h` | Complete |
| 7 | Graph Metrics (`graph_metrics_t`) | `include/subcubic.h` | Complete |
| 8 | Min-Plus Product (`min_plus_product_t`) | `include/subcubic.h` | Complete |
| 9 | OV Instance (`ov_instance_t`) | `include/subquadratic.h` | Complete |
| 10 | Edit Distance Instance (`edit_distance_instance_t`) | `include/subquadratic.h` | Complete |
| 11 | LCS Instance (`lcs_instance_t`) | `include/subquadratic.h` | Complete |
| 12 | Frechet Instance (`frechet_instance_t`) | `include/subquadratic.h` | Complete |
| 13 | DTW Instance (`dtw_instance_t`) | `include/subquadratic.h` | Complete |
| 14 | 3SUM Instance (`threesum_instance_t`) | `include/threesum.h` | Complete |
| 15 | Collinearity Instance (`collinearity_instance_t`) | `include/threesum.h` | Complete |
| 16 | Reduction Graph Node/Edge (`fg_node_t`, `fg_edge_t`) | `include/fine_grained_reduction.h` | Complete |

## L2: Core Concepts
| # | Concept | Implementation | Status |
|---|---------|---------------|--------|
| 1 | Subcubic Equivalence | `src/equiv_classes.c` | Complete |
| 2 | Subquadratic Equivalence | `src/equiv_classes.c` | Complete |
| 3 | 3SUM Equivalence | `src/equiv_classes.c` | Complete |
| 4 | Fine-Grained Reduction Types | `include/equiv_classes.h` | Complete |
| 5 | Equivalence Class Registration | `src/equiv_map.c` | Complete |
| 6 | Reduction Chain Computation | `src/equiv_classes.c` | Complete |
| 7 | Conjecture Verification | `src/equiv_classes.c` | Complete |
| 8 | Problem Classification | `src/equiv_classes.c` | Complete |

## L3: Mathematical Structures
| # | Structure | Implementation | Status |
|---|----------|---------------|--------|
| 1 | Distance Matrix (n x n doubles) | `src/subcubic.c` | Complete |
| 2 | Boolean Vector (bit-packed uint32_t) | `src/subquadratic.c` | Complete |
| 3 | DP Table (2D int32_t array) | `src/subquadratic.c` | Complete |
| 4 | Polygonal Curve (x,y arrays) | `src/subquadratic.c` | Complete |
| 5 | Integer Set (sorted int64_t array) | `src/threesum.c` | Complete |
| 6 | Reduction Graph (adjacency lists) | `src/fine_grained_reduction.c` | Complete |
| 7 | Min-Plus Semiring Operations | `src/subcubic.c` | Complete |

## L4: Fundamental Laws/Theorems
| # | Theorem | Implementation | Status |
|---|---------|---------------|--------|
| 1 | Subcubic Equivalence (Williams & Williams 2013) | `src/subcubic.c` | Complete |
| 2 | APSP <=> Min-Plus Product | `src/subcubic.c` | Complete |
| 3 | APSP <=> Negative Triangle | `src/subcubic.c` | Complete |
| 4 | SETH => OV Conjecture (Williams 2005) | `src/equiv_map.c` | Complete |
| 5 | OV => Edit Distance (Backurs & Indyk 2016) | `src/equiv_map.c` | Complete |
| 6 | OV => LCS/DTW (Abboud et al. 2015) | `src/equiv_map.c` | Complete |
| 7 | OV => Frechet (Bringmann 2014) | `src/equiv_map.c` | Complete |
| 8 | 3SUM-hardness of Collinearity | `src/threesum.c` | Complete |
| 9 | Transitivity of Fine-Grained Reductions | `src/fine_grained_reduction.c` | Complete |

## L5: Algorithms/Methods
| # | Algorithm | Time Complexity | Status |
|---|----------|---------|--------|
| 1 | Floyd-Warshall (APSP) | Theta(n^3) | Complete |
| 2 | Min-Plus Naive Product | Theta(n^3) | Complete |
| 3 | Min-Plus Fast (AGM 1997) | O(n^omega log M) | Complete (fallback) |
| 4 | Graph Metrics (Radius/Median/Diameter) | O(n^2) after APSP | Complete |
| 5 | OV Brute Force | O(n^2 d / wordsize) | Complete |
| 6 | OV Williams (2005) | O(n^{2-1/O(log d/log n)}) | Complete |
| 7 | Edit Distance DP | O(n*m) | Complete |
| 8 | Edit Distance Linear Space | O(n*m), O(min(n,m)) space | Complete |
| 9 | LCS DP | O(n*m) | Complete |
| 10 | LCS Hirschberg | O(n*m), O(min(n,m)) space | Complete |
| 11 | Frechet Distance | O(n*m) | Complete |
| 12 | DTW Standard | O(n*m) | Complete |
| 13 | DTW Sakoe-Chiba | O((n+m)*w) | Complete |
| 14 | 3SUM Sort+Two-Pointer | O(n^2) | Complete |
| 15 | Collinearity via 3SUM | O(n^2) | Complete |
| 16 | Tarjan SCC on Reduction Graph | O(V+E) | Complete |
| 17 | Floyd-Warshall Transitive Closure | O(V^3) | Complete |
| 18 | Dijkstra on Reduction Graph | O(V^2) | Complete |

## L6: Canonical Problems
| # | Problem | Classification | Status |
|---|---------|---------------|--------|
| 1 | APSP | Subcubic Canonical | Complete |
| 2 | Negative Triangle | Subcubic-Equivalent | Complete |
| 3 | Min-Plus Product | Subcubic-Equivalent | Complete |
| 4 | Graph Radius/Median/Diameter | Subcubic-Equivalent | Complete |
| 5 | Orthogonal Vectors (OV) | Subquadratic Canonical | Complete |
| 6 | Edit Distance | Subquadratic-Equivalent | Complete |
| 7 | LCS | Subquadratic-Equivalent | Complete |
| 8 | Frechet Distance | Subquadratic-Equivalent | Complete |
| 9 | DTW | Subquadratic-Equivalent | Complete |
| 10 | 3SUM | 3SUM Canonical | Complete |
| 11 | Collinearity | 3SUM-Hard | Complete |

## L7: Applications (+3)
| # | Application | Implementation | Status |
|---|-----------|---------------|--------|
| 1 | DNA Sequence Alignment | `src/subquadratic.c` | Complete |
| 2 | Time Series Similarity Search (DTW) | `src/subquadratic.c` | Complete |
| 3 | Mars Rover Obstacle Avoidance (NASA) | `src/threesum.c` | Complete |
| 4 | Network Routing Analysis (APSP) | `examples/example_subcubic.c` | Complete |
| 5 | Computational Geometry (Collinearity) | `examples/example_threesum.c` | Complete |

## L8: Advanced Topics
| # | Topic | Implementation | Status |
|---|-------|---------------|--------|
| 1 | Cross-Class Relations | `src/equiv_classes.c` | Complete |
| 2 | Fine-Grained Landscape Analysis | `src/fine_grained_reduction.c` | Complete |
| 3 | Barrier-Breaking Detection | `src/fine_grained_reduction.c` | Complete |
| 4 | Reduction Composition | `src/fine_grained_reduction.c` | Complete |

## L9: Research Frontiers
| # | Topic | Reference | Status |
|---|-------|-----------|--------|
| 1 | Subcubic Conjecture (Is APSP in O(n^{3-eps})?) | `src/equiv_map.c` | Partial (documented) |
| 2 | OV Conjecture and SETH | `src/equiv_map.c` | Partial (documented) |
| 3 | 3SUM Conjecture | `src/equiv_map.c` | Partial (documented) |
| 4 | Combinatorial BMM Conjecture | `src/equiv_map.c` | Partial (documented) |
| 5 | Fine-Grained Approximability | `src/equiv_map.c` | Partial (documented) |