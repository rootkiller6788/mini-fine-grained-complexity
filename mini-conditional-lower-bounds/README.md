# mini-conditional-lower-bounds

## Module Status: COMPLETE ✅

Methodology of conditional lower bounds in fine-grained complexity.

- **L1-L6**: Complete (all core definitions, concepts, structures, theorems, algorithms, and canonical problems)
- **L7**: Partial+ (3 applications: dynamic string matching, reachability oracle, regex matching)
- **L8**: Partial+ (2/4 advanced topics: subcubic equivalence classes, NSETH)
- **L9**: Partial (documented: quantum fine-grained, meta-complexity, barriers)

## Line Count
include/ + src/ = 3713 lines (minimum: 3000)

## Core Definitions (L1)
- `Hypothesis`: Formal representation of a fine-grained hardness conjecture
- `ConditionalLowerBound`: "If H then P requires time T(n)" statement
- `HypothesisDatabase`: Collection of hypotheses and their implications
- `ReductionWeb`: Directed graph of fine-grained reductions
- `EquivalenceClass`: Set of mutually reducible problems
- `FineGrainedTimeBound`: Time complexity with exponent and polylog
- `SatInstance`: CNF formula representation for SETH
- `ThreeSumInstance`, `ApspInstance`: Canonical problem instances
- `VectorSet`: Binary vectors for OV problem

## Core Theorems (L4)
- **SETH → OV** (Williams, 2005): If OV is in O(n^{2-eps}), then SETH is false
- **SETH → Edit Distance** (Backurs-Indyk, 2015): Edit Distance requires n^{2-o(1)}
- **SETH → LCS** (Abboud-Backurs-Williams, 2015): LCS requires n^{2-o(1)}
- **SETH → Diameter** (Roditty-Williams, 2013): Diameter requires n^{2-o(1)}
- **APSP subcubic equivalence** (Williams-Williams, 2010): 32+ problems require n^{3-o(1)}
- **3SUM quadratic lower bound** (Patrascu, 2010): 3SUM requires n^{2-o(1)}

## Core Algorithms (L5)
- 4 SAT solvers: Brute-force, DPLL, Schoning (O(1.333^n)), PPSZ (O(1.307^n))
- 4 OV algorithms: Brute-force, Split-and-list, Light-bulb, FFT-based
- 2 3SUM algorithms: Brute-force (O(n^3)), Quadratic (O(n^2))
- 2 APSP algorithms: Floyd-Warshall, Dijkstra-all
- 4 Edit Distance: Wagner-Fischer, Needleman-Wunsch, Four-Russians, Ukkonen
- 3 LCS: Standard DP, Extract, Hunt-Szymanski
- Frechet distance, DTW distance
- Kosaraju SCC for equivalence class computation
- BFS-based lower bound chain construction

## Canonical Problems (L6)
- Orthogonal Vectors (OV), k-SAT, 3SUM, APSP
- Edit Distance, LCS, Frechet Distance, DTW
- Graph Diameter, Negative Triangle, Graph Radius, Betweenness Centrality
- Minimum Area Triangle, 3-Points-on-a-Line, Hitting Set, Subgraph Isomorphism
- Regular Expression Matching

## Reduction Web
```
SETH —→ OV —→ Edit Distance
  |      |——→ LCS
  |      |——→ Diameter
  |      |——→ Hitting Set
  |      |——→ Subgraph Isomorphism
  |      |——→ Regex Matching
  |——→ Frechet Distance
  |——→ DTW

3SUM Conjecture —→ 3-Points-on-a-Line
                —→ Minimum Area Triangle
                —→ Polygon Containment

APSP Conjecture —→ Negative Triangle
                —→ Graph Radius
                —→ Median
                —→ Betweenness Centrality
                —→ Replacement Paths
```

## Build
```bash
make clean && make        # builds and runs tests (C99, libc+libm)
make test                 # run test suite (39 tests)
```

## Files
```
include/
  condlb.h                  Core CLB data structures (154 lines)
  seth_bound.h              SETH definitions and SAT API (195 lines)
  ov_reduction.h            OV definitions and reductions (155 lines)
  reduction_web.h           Reduction web and equivalence classes (214 lines)
  threesum_apsp.h           3SUM and APSP definitions (184 lines)
  string_hardness.h         String problem definitions (181 lines)
src/
  condlb.c                  Core implementation (758 lines)
  seth_bound.c              SETH, SAT solvers (580+ lines)
  ov_reduction.c            OV algorithms, reductions (223 lines)
  reduction_web.c           Reduction web, SCC, equivalence (426 lines)
  threesum_apsp.c           3SUM, APSP implementations (389 lines)
  string_hardness.c         Edit distance, LCS, Frechet, DTW (368 lines)
tests/
  test_main.c               Comprehensive test suite (39 tests)
docs/
  knowledge-graph.md        L1-L9 knowledge coverage
  coverage-report.md        Completeness evaluation
  gap-report.md             Missing knowledge points
  course-alignment.md       Nine-university curriculum mapping
  course-tree.md            Prerequisite dependency tree
```

## References
- Williams (2018): *Fine-Grained Complexity* (textbook draft)
- Williams (2015): *Hardness of Easy Problems* (STOC 2015)
- Williams & Williams (2010): *Subcubic Equivalences* (FOCS 2010)
- Backurs & Indyk (2015): *Edit Distance Cannot Be Computed in Strongly Subquadratic Time* (STOC 2015)
- Abboud, Backurs, Williams (2015): *Tight Hardness Results for LCS* (FOCS 2015)
- Bringmann (2014): *Why Walking the Dog Takes Time* (FOCS 2014)
- Patrascu (2010): *Towards Polynomial Lower Bounds for Dynamic Problems* (STOC 2010)
- Impagliazzo & Paturi (2001): *On the Complexity of k-SAT* (JCSS 2001)
- MIT 6.841, Stanford CS358/CS354, Berkeley CS278
