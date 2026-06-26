# mini-equivalence-classes

**Fine-Grained Equivalence Classes of Problems**

## Module Status: COMPLETE

| Level | Status |
|-------|--------|
| L1 Definitions | Complete (16 items) |
| L2 Core Concepts | Complete (8 items) |
| L3 Math Structures | Complete (7 items) |
| L4 Fundamental Laws | Complete (9 items) |
| L5 Algorithms | Complete (18 items) |
| L6 Canonical Problems | Complete (11 items) |
| L7 Applications | Complete (5 items) |
| L8 Advanced Topics | Complete (4 items) |
| L9 Research Frontiers | Partial (5 documented) |

---

## Core Definitions (L1)

- `equiv_class_t`: Equivalence class (canonical problem, threshold exponent)
- `problem_descriptor_t`: Problem metadata (name, class, best known exponent)
- `reduction_descriptor_t`: Fine-grained reduction between two problems
- `fg_reduction_graph_t`: Directed graph of problems and reductions
- `apsp_instance_t`: All-Pairs Shortest Paths instance (distance matrix)
- `neg_triangle_instance_t`: Negative Triangle detection instance
- `graph_metrics_t`: Graph radius, diameter, median, eccentricity
- `min_plus_product_t`: Min-Plus matrix multiplication
- `ov_instance_t`: Orthogonal Vectors (bit-packed vectors)
- `edit_distance_instance_t`: Edit Distance (Levenshtein) instance
- `lcs_instance_t`: Longest Common Subsequence instance
- `frechet_instance_t`: Frechet Distance between polygonal curves
- `dtw_instance_t`: Dynamic Time Warping instance
- `threesum_instance_t`: 3SUM instance
- `collinearity_instance_t`: Collinearity detection (3SUM-hard)

## Core Theorems (L4)

- **Subcubic Equivalence** (Williams & Williams 2013, STOC):
  APSP, Negative Triangle, Min-Plus, Radius, Median, Diameter are
  subcubic-equivalent under O(n^{3-eps}) reductions

- **SETH => OV Conjecture** (Williams 2005, ICALP):
  SETH => OV requires n^{2-o(1)} time when d = omega(log n)

- **OV => Edit Distance** (Backurs & Indyk 2016, STOC):
  Edit Distance requires n^{2-o(1)} time under SETH

- **OV => LCS, DTW** (Abboud, Backurs, Williams 2015, FOCS):
  LCS and DTW require n^{2-o(1)} time under SETH

- **OV => Frechet** (Bringmann 2014, FOCS):
  Frechet Distance requires n^{2-o(1)} time under SETH

- **3SUM Conjecture** (Gajentaan & Overmars 1995):
  3SUM requires n^{2-o(1)} time (weaker, independent of SETH)

## Equivalence Classes

| Class | Canonical | Threshold | Conjectured LB | Members |
|-------|-----------|-----------|---------------|---------|
| Subcubic | APSP | n^3 | n^{3-o(1)} | 8+ |
| Subquadratic | OV | n^2 | n^{2-o(1)} (SETH) | 7+ |
| 3SUM | 3SUM | n^2 | n^{2-o(1)} | 3+ |

## Key Algorithms (L5)

| Algorithm | Time | Class |
|-----------|------|-------|
| Floyd-Warshall (APSP) | Theta(n^3) | Subcubic |
| Min-Plus Naive | Theta(n^3) | Subcubic |
| Graph Metrics | O(n^2) after APSP | Subcubic |
| OV Brute Force | O(n^2 d/w) | Subquadratic |
| Edit Distance DP | O(nm) | Subquadratic |
| LCS DP / Hirschberg | O(nm) | Subquadratic |
| Frechet Distance | O(nm) | Subquadratic |
| DTW Standard | O(nm) | Subquadratic |
| 3SUM Sort+Two-Pointer | O(n^2) | 3SUM |

## Applications (L7)

1. DNA Sequence Alignment (Edit Distance)
2. Time Series Search (DTW)
3. Mars Rover Obstacle Planning (NASA) (Collinearity)
4. Network Routing (APSP)
5. Computational Geometry (3SUM-hard problems)

## Course Mapping

- MIT 6.890: Fine-Grained Complexity
- Stanford CS254: Lower Bounds
- CMU 15-751: Theorist's Toolkit
- Berkeley CS278: Computational Complexity
- Princeton COS 522: Complexity Theory
- Caltech CS 154: Limits of Computation
- Cambridge Part III: Advanced Complexity
- Oxford: Advanced Complexity Theory
- ETH 263-4650: Advanced Complexity

## Build & Test

```
make            # build library, examples, benchmarks
make test       # run all unit tests (40+ tests)
make clean      # clean build artifacts
make run-examples  # run all examples
make run-benches   # run benchmarks
```

## References

- Williams & Williams (2013), STOC
- Williams (2015), ICDT
- Backurs & Indyk (2016), STOC 2015
- Abboud, Backurs, Williams (2015), FOCS
- Bringmann (2014), FOCS
- Gajentaan & Overmars (1995), Comp. Geom.
- Gronlund & Pettie (2014), FOCS
- Williams (2018), SIGACT News