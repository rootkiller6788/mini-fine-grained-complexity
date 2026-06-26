# Knowledge Graph: mini-k-clique-hardness

## Nine-Level Knowledge Coverage

### L1: Definitions (Complete)
- k-Clique decision problem: Given G=(V,E), integer k, does G contain K_k?
- k-Clique search problem: Find a k-clique if one exists
- k-Clique counting: Count all k-cliques in G (#W[1]-complete)
- Max-Clique: Find a clique of maximum cardinality (NP-hard)
- W[1]-hardness definition (parameterized complexity class)
- FPT (Fixed-Parameter Tractability): solvable in f(k)*n^{O(1)} time
- Exponential Time Hypothesis (ETH): 3SAT requires 2^{Omega(n)} time
- Strong Exponential Time Hypothesis (SETH): lim_{k->inf} s_k = 1
- Graph complement, adjacency matrix, degree, neighbors

### L2: Core Concepts (Complete)
- Parameterized complexity: FPT vs XP vs W[1] hierarchy
- Conditional lower bounds: "If X, then Y requires Omega(n^{alpha})"
- Reduction-preserving hardness (FPT reductions)
- The trivial O(n^k) algorithm and its optimality under SETH
- Sub-exponential time and ETH violations
- Weft of circuits: W[t] hierarchy definition
- Parameter-preserving vs parameter-increasing reductions

### L3: Mathematical Structures (Complete)
- Graph as adjacency matrix (O(1) edge query, O(n^2) space)
- Combinatorial enumeration: C(n,k) k-subsets
- Ramsey theory: R(s,t) bounds, existence guarantees
- Turan's theorem: extremal graph theory bounds
- Degeneracy ordering and core decomposition
- Boolean circuits with bounded fan-in/weft
- Color-coding: hash families, k-perfect hash functions
- Nesetril-Poljak tensor product graph construction

### L4: Fundamental Laws / Theorems (Complete)
- W[1]-completeness of k-Clique (Downey & Fellows 1995)
- ETH => k-Clique requires n^{Omega(k)} time (Chen et al. 2006)
- SETH => k-Clique requires n^{k-o(1)} time (Lokshtanov et al. 2011)
- Monien's theorem: k-Clique in degree-d graphs is FPT
- Moon-Moser: max 3^{n/3} maximal cliques in an n-vertex graph
- Turan's theorem: ex(n, K_{r+1}) = (1-1/r) * n^2/2
- k-Clique no polynomial kernel unless NP subseteq coNP/poly
- Sparsification Lemma (Impagliazzo-Paturi-Zane 2001)
- Nesetril-Poljak: k-Clique reduces to triangle detection in tensor graph

### L5: Algorithms / Methods (Complete)
- Brute-force enumeration: O(n^k) by checking all subsets
- Color-coding (Alon-Yuster-Zwick 1995): O((2e)^k * n^2)
- Derandomized color-coding via k-perfect hash families
- Bron-Kerbosch with pivot: O(3^{n/3}) for maximal cliques
- Branch-and-bound with coloring bound (Tomita et al.)
- Matrix multiplication: O(n^{omega * ceil(k/3)})
- Degeneracy-based clique counting (Chiba-Nishizeki 1985)
- Greedy graph coloring for upper bounds
- Inclusion-exclusion counting

### L6: Canonical Problems (Complete)
- k-Clique (Decision, Search, Counting, Enumeration)
- k-Independent Set (via complement graph)
- k-Vertex Cover (via complement + size transformation)
- 3SAT -> k-Clique FPT reduction (W[1]-hardness proof)
- Maximum Clique (global optimization variant)
- Triangle Detection (k=3 special case)
- Parameterized Clique in bounded-degree/planar/bipartite graphs

### L7: Applications (Partial+)
- Social network community detection (clique = fully connected community)
- Bioinformatics motif finding (detecting dense interaction subgraphs)
- Compiler optimization (register allocation via graph coloring)
- Random graph models (Erdos-Renyi G(n,p)) for benchmarking
- SAT solver heuristics (clique in conflict graphs)

### L8: Advanced Topics (Partial+)
- Fine-grained equivalence: k-Clique and matrix multiplication
- Gap-ETH consequences for approximate clique
- W-hierarchy: weft analysis, circuit depth and parameterized hardness
- Kernelization lower bounds: no polynomial kernel for k-Clique
- Bounded search tree analysis for FPT algorithms
- Derandomization via explicit hash families (Naor-Schulman-Srinivasan)
- Tensor power graph constructions (Nesetril-Poljak)

### L9: Research Frontiers (Partial+)
- Is W[1] = FPT? (major open problem)
- Gap-ETH and hardness of approximation
- Average-case complexity of k-Clique
- Planted clique problem: finding hidden cliques in random graphs
- Quantum algorithms for triangle/k-clique detection
- Sub-exponential algorithms for special graph classes
