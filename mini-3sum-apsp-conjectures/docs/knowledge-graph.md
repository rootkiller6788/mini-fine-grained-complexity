# Knowledge Graph: 3SUM & APSP Conjectures

## L1: Definitions
- 3SUM problem: Given n integers, find a+b+c=0
- APSP problem: Shortest paths between all vertex pairs
- Min-plus matrix product: C[i][j]=min_k{A[i][k]+B[k][j]}
- Fine-grained reduction: preserves exact polynomial exponent
- 3SUM-hard: Problem requiring Omega(n^2) assuming 3SUM conjecture
- APSP-hard: Problem requiring Omega(n^3) assuming APSP conjecture

## L2: Core Concepts
- Quadratic barrier (3SUM): O(n^{2-epsilon}) conjectured impossible
- Cubic barrier (APSP): O(n^{3-epsilon}) conjectured impossible
- Subcubic equivalence: 3SUM<->APSP under subcubic reductions
- Conditional lower bounds: If conjecture holds, problem is hard
- Hardness classification: 3SUM-hard, APSP-hard, SETH-hard families

## L3: Mathematical Structures
- Integer sets and triples (3SUM)
- Weighted directed graphs (APSP)
- Distance matrix and adjacency matrix
- Min-plus (tropical) semiring: (R U {inf}, min, +)
- Hash table with open addressing (3SUM optimization)
- Algebraic decision tree model (3SUM lower bounds)

## L4: Fundamental Theorems
- 3SUM Conjecture (Gajentaan-Overmars 1995)
- APSP Conjecture (Vassilevska-Williams 2009)
- Subcubic Equivalence Theorem (VW 2009): 3SUM and APSP are equivalent
- Erickson Lower Bound (1995): Omega(n^2) in 3-linear decision tree
- Gronlund-Pettie Theorem (2014): O(n^2(log log n/log n)^{2/3})
- Gold-Sharir Theorem (2017): O(n^2 log log n / log n)

## L5: Algorithms
- Naive 3SUM: O(n^3) triple loop
- Hash-based 3SUM: O(n^2) expected time
- Sort+bsearch 3SUM: O(n^2 log n) deterministic
- Two-pointer 3SUM: O(n^2) deterministic
- Gronlund-Pettie 2014: O(n^2/(log n)^{2/3}) breakthrough
- Floyd-Warshall APSP: O(n^3) dynamic programming
- Johnson APSP: O(nm + n^2 log n) for sparse graphs
- Dijkstra-all APSP: O(n(m + n log n)) non-negative weights
- Bellman-Ford-all APSP: O(n^2 m) with negative edges

## L6: Canonical Problems
- 3SUM (baseline subquadratic problem)
- Collinearity Detection (3SUM-hard)
- Minimum Area Triangle (3SUM-hard)
- APSP (baseline subcubic problem)
- Min-Plus Matrix Product (APSP-equivalent)
- Negative Triangle Detection (APSP-equivalent)
- Graph Radius (APSP-equivalent)
- Graph Median (APSP-equivalent)
- Betweenness Centrality (APSP-equivalent)

## L7: Applications
- Computational geometry: collinearity testing
- Network routing: shortest-path computation
- Traffic analysis: bottleneck detection
- Social network analysis: betweenness centrality
- Database queries: join optimization (3SUM related)
- Motion planning: configuration space obstacles

## L8: Advanced Topics
- Subcubic equivalence hexagon (Williams & Williams 2010)
- Polynomial method for 3SUM (Gajentaan-Overmars framework)
- Fine-grained complexity classes beyond 3SUM/APSP
- Real-RAM vs integer RAM model distinctions
- Word-RAM parallelism (Gronlund-Pettie technique)
- Deterministic hashing (Fredman-Komlos-Szemeredi)

## L9: Research Frontiers
- Quantum algorithms for 3SUM and APSP
- Beating n^2 for 3SUM in the decision tree model
- APSP in O(n^{2.999}) - still open
- Min-plus product faster than O(n^3) - still open
- Unifying SETH, 3SUM, and APSP hardness
- Meta-complexity of fine-grained reductions
