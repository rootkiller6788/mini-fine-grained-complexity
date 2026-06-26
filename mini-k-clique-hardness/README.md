# mini-k-clique-hardness

## Status: COMPLETE ✅

k-Clique hardness under ETH/SETH. W[1]-completeness and conditional lower bounds.

**L1-L4: Complete | L5-L6: Complete | L7-L9: Partial+**

### Core Definitions
- **k-Clique Decision**: Does G contain a k-clique? (W[1]-complete)
- **k-Clique Search**: Find a k-clique if one exists
- **k-Clique Counting**: Count all k-cliques (#W[1]-complete)
- **Max-Clique**: Find the maximum clique (NP-hard)
- **ETH**: Exponential Time Hypothesis (3SAT requires 2^{Ω(n)} time)
- **SETH**: Strong ETH (kSAT requires 2^{(1-o(1))n} time)
- **FPT**: Fixed-Parameter Tractability (f(k)·n^{O(1)} time)
- **W[1]**: First level of parameterized intractability hierarchy

### Core Theorems
- **W[1]-completeness** (Downey & Fellows 1995): k-Clique is W[1]-complete
- **ETH ⇒ n^{Ω(k)}** (Chen et al. 2006): k-Clique requires n^{Ω(k)} time under ETH
- **SETH ⇒ n^{k-o(1)}** (Lokshtanov et al. 2011): k-Clique requires n^{k-o(1)} under SETH
- **Monien's Theorem** (1985): k-Clique in degree-d graphs is FPT: O(k²·d^k·n)
- **Moon-Moser** (1965): At most 3^{n/3} maximal cliques in n-vertex graph
- **Turán's Theorem** (1941): ex(n, K_{r+1}) = (1-1/r)·n²/2
- **Nesetril-Poljak** (1985): k-Clique reduces to triangle in tensor graph

### Core Algorithms
- **Brute Force**: O(n^k) — trivial, optimal under SETH
- **Color-Coding** (Alon-Yuster-Zwick 1995): O((2e)^k · n²) — FPT
- **Bron-Kerbosch** (1973): O(3^{n/3}) — optimal maximal clique enumeration
- **Branch & Bound** (Tomita et al. 2003): O(1.442^n) worst-case
- **Matrix Multiplication**: O(n^{ω·⌈k/3⌉}) — beats brute force for small k
- **Degeneracy Counting** (Chiba-Nishizeki 1985): O(k·d^{k-2}·m)

### Canonical Problems
- k-Clique (Decision, Search, Counting, Enumeration)
- k-Independent Set (via complement)
- k-Vertex Cover (via complement + size)
- 3SAT → k-Clique (W[1]-hardness reduction)

### Nine-School Course Mapping
- MIT 6.841: Fine-grained complexity, SETH
- Stanford CS254: Conditional lower bounds
- Berkeley CS278: W-hierarchy, parameterized complexity
- CMU 15-855: Fine-grained equivalence
- Princeton COS 522: ETH/SETH
- Caltech CS 154: Fine-grained analysis
- Cambridge Part III: Parameterized complexity
- Oxford: FPT algorithms
- ETH 263-4650: SETH lower bounds

### Nine-Level Knowledge Coverage

| Level | Status | Core Topics |
|-------|--------|------------|
| L1 Definitions | Complete | 11 struct types, graph/clique/ETH/SETH definitions |
| L2 Core Concepts | Complete | FPT, W[1], conditional lower bounds |
| L3 Math Structures | Complete | Adjacency matrix, Turan bounds, Ramsey numbers |
| L4 Fundamental Laws | Complete | W[1]-completeness, ETH/SETH lower bounds |
| L5 Algorithms | Complete | 6 major algorithm families implemented |
| L6 Canonical Problems | Complete | k-Clique, k-IS, k-VC, 3SAT |
| L7 Applications | Partial+ | Social networks, bioinformatics, compiler |
| L8 Advanced Topics | Partial+ | Gap-ETH, kernel bounds, derandomization |
| L9 Research Frontiers | Partial+ | W[1]=FPT, quantum, planted clique |

### Build
```
make          # Build library, tests, examples, benchmarks
make test     # Run test suite
make run-examples  # Run all examples
```

### Files
- `include/`: 6 header files (~1450 lines)
- `src/`: 9 C files + 1 Lean file (~3200 lines)
- `tests/`: Comprehensive test suite
- `examples/`: 4 executable examples
- `docs/`: 5 knowledge documents
- `benches/`: Performance benchmarks

### References
- Downey & Fellows, "Parameterized Complexity" (1999)
- Chen, Huang, Kanj, Xia, "Tight Lower Bounds" (JCSS 2006)
- Lokshtanov, Marx, Saurabh, "Lower Bounds for k-Clique" (SODA 2011)
- Alon, Yuster, Zwick, "Color-Coding" (JACM 1995)
- Impagliazzo & Paturi, "Complexity of k-SAT" (JCSS 2001)
