# Course Tree: mini-k-clique-hardness

## Prerequisite Knowledge Dependencies

```
                            [Computability Basics]
                                     |
                            [P vs NP, NP-Completeness]
                                     |
                    +----------------+----------------+
                    |                |                |
            [Graph Theory]    [3SAT, kSAT]    [Combinatorics]
                    |                |                |
                    +--------+-------+--------+-------+
                             |                |
                     [Parameterized      [Fine-Grained
                      Complexity]         Complexity]
                             |                |
                             +-------+--------+
                                     |
                          [k-Clique Hardness]
                                     |
                    +--------+-------+-------+--------+
                    |        |               |        |
              [ETH/SETH]  [W-Hierarchy]  [Color-   [Matrix
                           [Kernel LB]    Coding]   Mult]
```

## Dependency Details

### Direct Prerequisites
1. **Graph Theory** (L3): Adjacency, cliques, independent sets, complements
2. **3SAT/kSAT** (L6): Canonical NP-complete problems, clause structure
3. **Combinatorics** (L3): Binomial coefficients, subset enumeration, Ramsey
4. **Parameterized Complexity** (L2): FPT, W-hierarchy, kernelization
5. **Fine-Grained Complexity** (L2): ETH, SETH, conditional lower bounds

### Indirect Prerequisites
- **Complexity Theory** (L1-L2): P, NP, polynomial-time reductions
- **Algorithm Design** (L5): Dynamic programming, search algorithms
- **Linear Algebra** (L3): Matrix multiplication, adjacency matrices

## Knowledge Flow

```
L1: Definitions ? structs, types (graph_t, clique_t, eth_params_t)
  ?
L2: Core Concepts ? FPT framework, conditional lower bound encoding
  ?
L3: Math Structures ? Graph operations, degeneracy, core decomposition
  ?
L4: Fundamental Laws ? ETH?n^{?(k)}, SETH?n^{k-o(1)}, W[1]-completeness
  ?
L5: Algorithms ? Brute-force, BK, color-coding, B&B, matrix multiplication
  ?
L6: Canonical Problems ? k-Clique, k-IS, k-VC, 3SAT reduction
  ?
L7: Applications ? Social networks, bioinformatics, compiler optimization
  ?
L8: Advanced ? Gap-ETH, kernel lower bounds, derandomization
  ?
L9: Frontiers ? W[1]=FPT, quantum algorithms, planted clique
```

## Module Dependencies

This module depends on:
- `mini-complexity-foundations/mini-p-np-np-completeness` (NP, reductions)
- `mini-complexity-foundations/mini-reductions-completeness` (SAT, 3SAT)
- `mini-circuit-complexity/mini-boolean-circuits-model` (circuit weft)

This module is a prerequisite for:
- `mini-fine-grained-complexity/mini-conditional-lower-bounds` (methodology)
- `mini-fine-grained-complexity/mini-equivalence-classes` (SETH-based classes)
- `mini-fine-grained-complexity/mini-polynomial-method-consequences`
