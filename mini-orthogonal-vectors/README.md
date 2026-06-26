# mini-orthogonal-vectors

## Status: COMPLETE ✅

Orthogonal Vectors problem and its connection to SETH.

Williams (2005): SETH → OV conjecture. Abboud-Williams-Yu (2015): Equivalence class.

## Module Stats
- **include/ + src/ lines**: 3044 (≥ 3000 ✓)
- **Tests passed**: 84/84
- **Compilation**: `make` passes with -Wall -Wextra
- **Headers**: 2 (ov.h, ov_reductions.h)
- **Source files**: 10 (.c files in src/)
- **Docs**: 5 knowledge documents in docs/

## Nine-Layer Knowledge Coverage

| Level | Status    | Key Content |
|-------|-----------|-------------|
| L1    | Complete  | binary_vector_t, vector_set_t, ov_instance_t, 20+ operations |
| L2    | Complete  | OV Conjecture formalization, 3-regime classification |
| L3    | Complete  | Packed bit ops, GF(2) algebra, Gaussian elimination |
| L4    | Complete  | Williams' Theorem (SETH→OV), lower bounds, Sparsification |
| L5    | Complete  | 8 algorithms: Brute, Williams, MITM, 4-Russians, MM, Poly, LSH |
| L6    | Complete  | 5 reductions: PM, Diameter, EditDist, LCS, SubsetSum |
| L7    | Complete  | 4 apps: Database, CRISPR, CodeSim, CollaborativeFilter |
| L8    | Partial+  | Dynamic OV, Approx OV, CommComplexity, Difficulty estimation |
| L9    | Partial   | 3SUM implication, non-binary OV, equivalence class |

## Core Definitions (L1)
- `binary_vector_t` — Binary vector in {0,1}^d, packed bits (32 bits/word)
- `vector_set_t` — Collection of n vectors of dimension d
- `ov_instance_t` — Pair of vector sets (A, B)
- `ov_result_t` — Query result with indices, ops count, timing
- Orthogonality: <a,b> = Σ a_i·b_i = 0

## Core Theorems (L4)
- **Williams' Theorem (2005)**: SETH ⇒ OV Conjecture
  - k-SAT requiring 2^{(1-ε)n} for all k ⇒ OV requires n^{2-o(1)} for d=ω(log n)
- **OV Lower Bound**: d = ω(log n) ⇒ Ω(n^{2-o(1)})
- **OV → PM/ED/LCS/Diameter/SubsetSum**: Subquadratic reductions
- **Polynomial Method**: Degree-d polynomial detection requires n^{2-o(1)} for d=ω(log n)

## Core Algorithms (L5)
| Algorithm | Time | Space | When Best |
|-----------|------|-------|-----------|
| Brute Force | O(n²·d) | O(1) | Always works |
| Williams | n^{2-1/O(c)} | O(n·d) | d = c·log n |
| Meet-in-Middle | O(n·2^{d/2}) | O(n+2^{d/2}) | d ≤ 24 |
| Four-Russians | O(n²·d/log n) | O(2^b) | d = Θ(log n) |
| Matrix Mult | O(n^ω) | O(n²) | ω<2.373 |
| Polynomial | n^{2-ε} | O(n·d) | d = c·log n |
| Light-Sparse | O(n²·δ·d) | O(n·d) | δ=o(1) |
| LSH | n^{1+ρ} | O(n·d) | Approximation |

## Canonical Problems (L6)
OV is subquadratic-equivalent to:
- Pattern Matching with Wildcards
- Edit Distance (Levenshtein)
- Longest Common Subsequence
- Frechet Distance
- Dynamic Time Warping
- Graph Diameter (sparse)
- Hamming Nearest Neighbor
- Subset Sum (bounded integers)

## Nine-School Course Mapping
MIT 6.841, Stanford CS254, Berkeley CS278, CMU 15-855,
Princeton COS 551, Caltech CS 154, Cambridge Part III,
Oxford Advanced Complexity, ETH 263-4650.

Reference: Arora & Barak (2009), Abboud-Williams-Yu (2015), Williams (2005).

## Build & Test
```
make          # build library + tests
make test     # run tests (84/84 pass)
make clean    # clean build
make examples # build examples
make benches  # build benchmarks
```

See mini-seth-strong-eth for SETH definitions.
