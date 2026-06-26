# doc-3: Exponential-Time SAT Algorithms

## L5: Algorithms and Their Complexity

### 3.1 Algorithmic Landscape for k-SAT

The quest for the fastest k-SAT algorithm directly determines our understanding of SETH. Each improvement to the exponent s_k pushes the SETH limit lower.

| Algorithm | Year | 3-SAT Exponent | k-SAT Time | Key Technique |
|-----------|------|----------------|------------|---------------|
| Brute Force | — | 1.0000 | O(2^n) | Enumeration |
| DPLL | 1962 | 1.0000 | O(2^n) | Backtracking + unit propagation |
| Monien-Speckenmeyer | 1980 | ~0.667 | O(1.6181^n) | Autarky, splitting |
| Schöning | 1999 | ~0.415 | O((2(k-1)/k)^n) | Random walk |
| Paturi-Pudlak-Zane | 1998 | ~0.386 | O(2^{(1-μ_k)n}) | Resolution + random order |
| PPSZ | 2005 | 0.3863 | O(2^{(1-1/O(k))n}) | Bounded-width resolution |

### 3.2 Brute Force (L5)

**Time:** O(2^n)
**Space:** O(n)
**Method:** Enumerate all 2^n assignments, check each against formula.
**Role:** Baseline for all SETH lower bounds. Any algorithm beating 2^n for general CNF would refute SETH.

### 3.3 Schöning's Random Walk Algorithm (1999)

**Core Idea:** Start from a random assignment. If unsatisfied, pick a random unsatisfied clause and flip a random literal in it. Repeat.

**Expected time for k-SAT:** O((2(k-1)/k)^n)

For 3-SAT: (2·2/3)^n = (4/3)^n = 2^{0.415n}

**Success probability per trial:** (k/(2(k-1)))^n

**Implementation:**
```c
sat_result_t sat_schoening(const cnf_formula_t *f,
                           int32_t max_tries, int32_t max_flips);
```

### 3.4 PPSZ Algorithm (2005) — Current Best

**Core Idea:**
1. Pick a random ordering of variables
2. For each variable in order:
   - If its value is forced by a small clause (bounded-width resolution), set it accordingly
   - Otherwise, guess randomly
3. Use bounded-width resolution preprocessing to identify forced assignments

**3-SAT upper bound:** O(2^{0.3863n})
**k-SAT upper bound:** O(2^{(1-μ_k)n}) where μ_k = Ω(1/k)

**Key insight:** Bounded-width resolution (width ≤ 3) catches most forced variables, dramatically reducing the search space.

### 3.5 CDCL: Conflict-Driven Clause Learning (1996-present)

**Backbone of modern SAT solvers (MiniSAT, Glucose, CaDiCaL):**

1. **Boolean Constraint Propagation (BCP):** Two-watched-literal scheme for efficient unit propagation
2. **Variable State Independent Decaying Sum (VSIDS):** Activity-based branching heuristic
3. **Conflict Analysis:** First UIP (Unique Implication Point) learning
4. **Learned Clauses:** Prevent repeating the same conflict
5. **Restarts:** Periodically reset to escape local minima
6. **Clause Deletion:** Remove less useful learned clauses

**Theoretical worst case:** O(2^n), but empirically solves millions of variables.
**SETH significance:** CDCL does NOT refute SETH — its worst case is still exponential.

### 3.6 Algorithm Comparison for 3-SAT (n=50)

| Solver | Exponent | Operations (approx) | Feasible? |
|--------|----------|---------------------|-----------|
| Brute Force | 1.0000 | 2^50 ≈ 10^15 | No |
| DPLL | 1.0000 | 2^50 ≈ 10^15 | No |
| Schöning | 0.4150 | 2^20.75 ≈ 2×10^6 | Yes |
| PPSZ | 0.3863 | 2^19.3 ≈ 6×10^5 | Yes |

## References

1. Schöning, U. (1999). "A Probabilistic Algorithm for k-SAT and Constraint Satisfaction Problems." FOCS 1999.
2. Paturi, R., Pudlak, P., Saks, M.E., Zane, F. (2005). JACM 52(3):337-364.
3. Davis, M., Putnam, H. (1960). "A Computing Procedure for Quantification Theory." JACM 7(3):201-215.
4. Davis, M., Logemann, G., Loveland, D. (1962). CACM 5(7):394-397.
5. Marques-Silva, J., Sakallah, K. (1996). "GRASP — A New Search Algorithm for Satisfiability." ICCAD 1996.
