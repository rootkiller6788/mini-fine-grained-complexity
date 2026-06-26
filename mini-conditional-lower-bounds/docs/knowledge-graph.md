# Knowledge Graph — mini-conditional-lower-bounds

## L1: Definitions
- Conditional Lower Bound (CLB): "If Hypothesis H, then Problem P requires time T(n)"
- Fine-Grained Reduction: reduction preserving sub-polynomial factors
- SETH (Strong Exponential Time Hypothesis): lim_{k->inf} s_k = 1
- ETH (Exponential Time Hypothesis): 3-SAT requires 2^{Omega(n)} time
- OVC (Orthogonal Vectors Conjecture): OV requires n^{2-o(1)} time
- 3SUM Conjecture: 3SUM requires n^{2-o(1)} time
- APSP Conjecture: APSP requires n^{3-o(1)} time
- NSETH (Nondeterministic SETH)
- Hypothesis Database: collection of hypotheses and their implications
- Reduction Web: directed graph of fine-grained reductions
- Equivalence Class: set of mutually reducible problems
- FineGrainedTimeBound: time complexity bound with exponents
- ReductionEdge: weighted directed edge in reduction web
- ProblemNode: computational problem in reduction web

## L2: Core Concepts
- Conditional lower bound methodology
- Fine-grained complexity classes (subcubic, subquadratic)
- Hypothesis-based hardness: basing hardness on unproven conjectures
- Reduction web structure and inference
- Tightness of lower bounds
- Polylog-factor preservation
- Randomized fine-grained reductions
- Hardness tier computation
- Subcubic equivalence class
- Subquadratic equivalence class
- 3SUM-hard problem class
- Lower bound chain construction

## L3: Mathematical Structures
- Hypothesis struct with implication DAG
- HypothesisDatabase with hypotheses and bounds
- ReductionWeb as adjacency-list graph
- ExtendedReductionWeb with problem metadata
- EquivalenceClass with member sets
- WeightedReduction with oracle calls and blowup
- LowerBoundChain: sequential reduction path
- FineGrainedTimeBound with exponent and polylog
- Kosaraju SCC algorithm for equivalence classes
- BFS-based path existence in reduction web
- Transitive closure computation

## L4: Fundamental Theorems
- SETH-to-OV reduction theorem (Williams, 2005)
- SETH-to-Edit Distance theorem (Backurs-Indyk, 2015)
- SETH-to-LCS theorem (Abboud-Backurs-Williams, 2015)
- SETH-to-Diameter theorem (Roditty-Williams, 2013)
- 3SUM quadratic lower bound theorem (Patrascu, 2010)
- APSP subcubic-equivalence theorem (Williams-Williams, 2010)
- Subcubic equivalence class size (32+ problems)
- SETH-to-Hitting Set (Abboud et al., 2016)
- Frechet Distance SETH hardness (Bringmann, 2014)
- DTW SETH hardness (Bringmann-Kunnemann, 2015)

## L5: Algorithms/Methods
- Brute-force OV detection: O(n^2 * d)
- Split-and-list OV: O~(n^{2-1/O(log c)}) (Williams, 2005)
- Light-bulb OV for sparse vectors: O(n^2 * s/d)
- FFT-based OV for small dimension: O(n * 2^{d/2})
- SAT brute-force: O(2^n * m)
- DPLL SAT solving with unit propagation and pure literals
- Schoning randomized k-SAT: O((2(k-1)/k)^n)
- PPSZ randomized k-SAT: O(1.307^n) for k=3
- O(n^2) 3SUM via sorting + two-pointer
- Floyd-Warshall APSP: O(n^3)
- Dijkstra-all APSP: O(n * (m + n log n))
- Wagner-Fischer edit distance DP: O(n*m)
- Four-Russians edit distance: O(n^2 / log^2 n)
- Ukkonen edit distance: O(n + d^2)
- LCS DP: O(n*m)
- Hunt-Szymanski LCS: O((R+n) log n)
- Frechet distance DP: O(n*m)
- DTW distance DP: O(n*m)
- Kosaraju SCC decomposition: O(V+E)
- Lower bound chain construction via BFS

## L6: Canonical Problems
- Orthogonal Vectors (OV)
- k-SAT (3-SAT, CNF-SAT)
- 3SUM
- All-Pairs Shortest Paths (APSP)
- Edit Distance (Levenshtein)
- Longest Common Subsequence (LCS)
- Frechet Distance
- Dynamic Time Warping (DTW)
- Graph Diameter
- Negative Triangle
- Graph Radius
- Betweenness Centrality
- Minimum Area Triangle
- 3-Points-on-a-Line (Collinearity)
- Hitting Set
- Subgraph Isomorphism
- Regular Expression Matching

## L7: Applications
- Dynamic string matching lower bounds (Abboud-Williams, 2014)
- Graph reachability oracle bounds
- Regular expression matching lower bounds (Backurs-Indyk, 2016)
- Computational geometry lower bounds via 3SUM
- Social network analysis lower bounds via APSP
- Bioinformatics sequence alignment lower bounds

## L8: Advanced Topics
- Subcubic equivalence class structure (Williams-Williams, 2018)
- Subquadratic equivalence class structure
- NSETH and its implications
- Polynomial method in fine-grained complexity
- Hardness in P via conditional lower bounds
- SETH-to-circuit-SAT reduction

## L9: Research Frontiers
- Fine-grained complexity beyond polynomial time
- Equivalence class completeness
- Quantum fine-grained complexity
- Meta-complexity and fine-grained hardness
- Barriers to proving unconditional lower bounds
