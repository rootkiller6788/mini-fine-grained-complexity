# mini-seth-strong-eth

**Strong Exponential Time Hypothesis (SETH) and Exponential Time Hypothesis (ETH)**

## Module Status: COMPLETE ✅

- L1-L6: Complete
- L7: Complete (3+ applications)
- L8: Complete (SETH implications, circuit lower bounds)
- L9: Partial (documented, empirical validation)

## Core Definitions (L1)

- `cnf_formula_t`: CNF formula representation (clauses, variables)
- `assignment_t`: Partial/full truth assignment
- `literal_t`: Variable with sign bit
- `clause_t`: Disjunction of literals

## Core Theorems (L4)

- **SETH**: ∀ε>0, ∃k s.t. k-SAT ∉ DTIME(2^{(1-ε)n})
- **ETH**: 3-SAT ∉ DTIME(2^{o(n)})
- **Sparsification Lemma** (IPZ 2001): Any k-CNF → OR of ≤2^{εn} k-CNFs with O(n) clauses
- **Isolation Lemma** (VV 1986): Unique-SAT reduction with prob ≥1/(4n)

## Core Algorithms (L5)

- Brute-force SAT: O(2^n) baseline
- DPLL (Davis-Putnam-Logemann-Loveland 1962): Backtracking + unit propagation
- Schöning (1999): Random walk, O((2(k-1)/k)^n) = O(1.334^n) for 3-SAT
- PPSZ (2005): Bounded-width resolution + random ordering, O(2^{0.3863n}) for 3-SAT
- CDCL: Conflict-driven clause learning (modern SAT backbone)

## Key Numbers

| k | s_k (SETH limit) | Best Algorithm |
|---|-----------------|----------------|
| 3 | 0.3863 | PPSZ: 2^{0.3863n} |
| 4 | 0.5548 | PPSZ: 2^{0.5548n} |
| 5 | 0.6500 | PPSZ |
| ∞ | 1.0000 | SETH conjectured |

## Course Mapping

- MIT 6.890: Fine-Grained Complexity
- Stanford CS254 §15: Lower Bounds
- CMU 15-751: Theorist's Toolkit

## Build & Test

```
make        # build library
make test   # run 13 tests (all passing)
make clean  # clean build
```
