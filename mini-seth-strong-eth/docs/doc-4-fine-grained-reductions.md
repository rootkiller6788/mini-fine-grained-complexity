# doc-4: Fine-Grained Reductions from SETH

## L6-L7: Conditional Lower Bounds

### 4.1 The Reduction Framework

A **fine-grained reduction** from problem A to problem B shows: if B can be solved in time T(n), then A can be solved in time T(n)^{O(1)}. The key is **preserving the exact exponent**, not just showing NP-hardness.

### 4.2 Orthogonal Vectors (OV) — The Foundational Reduction

**Problem:** Given two sets A, B of N binary vectors in {0,1}^d, are there a ∈ A, b ∈ B with ⟨a, b⟩ = 0?

**SETH → OV (Williams 2005):**
1. Split n variables into two groups of n/2
2. For each partial assignment to group 1, create a vector encoding which clauses are satisfied
3. For each partial assignment to group 2, create a vector encoding which clauses still need satisfaction
4. Formula is SAT iff there exists an orthogonal pair

**Result:** OV instance with N = 2^{n/2} vectors of dimension d = m.
An O(N^{2-ε}) OV algorithm would solve k-SAT in 2^{(1-ε/2)n} time, refuting SETH.

**Conditional lower bound:** Under SETH, OV requires N^{2-o(1)} time for d = ω(log N).

### 4.3 Edit Distance — Quadratic-Time Hardness

**Problem:** Compute the Levenshtein distance between two strings of length N.

**SETH → Edit Distance (Backurs-Indyk 2015):**
- Encode SAT formula structure into strings
- Edit distance is small iff the two partial assignments together satisfy the formula
- String length: N = O(2^{n/2} · poly(n))

**Conditional lower bound:** Under SETH, Edit Distance requires N^{2-o(1)} time.
This explains why no O(N^{2-ε}) algorithm has been found despite decades of effort.

### 4.4 Longest Common Subsequence (LCS)

**Problem:** Find the length of the longest subsequence common to two strings of length N.

**SETH → LCS (Abboud-Backurs-Williams 2015):**
- Gadget construction from clause structure
- An O(N^{2-ε}) LCS algorithm would refute SETH

**Current best:** O(N² / log² N) via Four-Russians technique (Masek-Paterson 1980)

### 4.5 k-Dominating Set — Tight Parameterized Lower Bound

**Problem:** In a graph with n vertices, find k vertices that dominate all others.

**SETH → k-Dominating Set (Patrascu-Williams 2010):**
- Group n variables into k groups of size n/k
- Create a set system where each group becomes a "choice" gadget
- Reduction produces graph with n' = O(k · 2^{n/k}) vertices

**Conditional lower bound:** Under SETH, k-Dominating Set requires n^{k-o(1)} time.
This is **tight** — matches the O(n^{k+1}) trivial algorithm up to the o(1) in the exponent.

### 4.6 The Reduction Graph

Our reduction graph tracks the web of fine-grained reductions:

| From | To | Blowup | Reference |
|------|----|--------|-----------|
| k-SAT | OV | n^(0.5) | Williams 2005 |
| k-SAT | Edit Distance | n^(0.5) | Backurs-Indyk 2015 |
| k-SAT | LCS | n^(0.5) | ABW 2015 |
| k-SAT | k-Dom Set | n^(0.333) | PW 2010 |
| OV | k-Clique | n^(1.0) | Folklore |

### 4.7 Other SETH-Based Lower Bounds

| Problem | Lower Bound under SETH | Best Known UB |
|---------|----------------------|---------------|
| Frechet Distance | n^{2-o(1)} | O(n² log n) |
| Dynamic Connectivity | n^{1-o(1)} per op | O(log n) amortized |
| 3SUM | n^{2-o(1)} | O(n² / (log n)²) |
| APSP | n^{3-o(1)} | O(n³ / 2^{Ω(√log n)}) |

## References

1. Williams, R. (2005). "A new algorithm for optimal 2-constraint satisfaction and its implications." TCS 348(2):357-365.
2. Backurs, A., Indyk, P. (2015). "Edit Distance Cannot Be Computed in Strongly Subquadratic Time (unless SETH is false)." STOC 2015.
3. Abboud, A., Backurs, A., Williams, V.V. (2015). "Tight Hardness Results for LCS and Other Sequence Similarity Measures." FOCS 2015.
4. Patrascu, M., Williams, R. (2010). "On the Possibility of Faster SAT Algorithms." SODA 2010.
