# Knowledge Graph - mini-polynomial-method-consequences

## L1: Core Definitions (Complete)

| Concept | C Definition | Description |
|---------|-------------|-------------|
| Polynomial | polynomial_t | Multivariate polynomial over GF(2) or reals |
| Monomial | monomial_t | Product of variables (bitmask) |
| Term | term_t | Coefficient times monomial |
| GF(2) element | gf2_t | 2-element field element |
| Algebraic circuit | algebraic_circuit_t | DAG computing polynomials |
| OV instance | ov_instance_t | Set of Boolean vectors |
| CNF formula | williams_cnf_t | Conjunctive Normal Form |
| Assignment | williams_assign_t | Variable truth assignment |

## L2: Core Concepts (Complete)

- Polynomial degree, multi-linear polynomials
- Truth table to polynomial conversion
- Fourier expansion via FWHT
- Subset convolution (zeta/mobius transforms)
- Sparse polynomial representation

## L3: Mathematical Structures (Complete)

- GF(2) arithmetic: XOR/AND
- GF(p) arithmetic: modular add/mul/pow/inv
- Polynomial ring over Boolean variables
- Boolean cube and Fourier basis
- Subset lattice and ranked posets

## L4: Fundamental Laws (Complete)

- Zhegalkin theorem: unique GF(2) multi-linear polynomial
- Fourier uniqueness: unique spectral representation
- Razborov approximation: AC0 approximated by low-degree polynomials
- Smolensky lower bound: MOD_q not in AC0[p] for distinct primes
- Williams correctness: polynomial-method SAT algorithm

## L5: Algorithms (Complete)

- GF(2) Mobius transform O(n*2^n)
- Fast Walsh-Hadamard transform O(n*2^n)
- Zeta/Mobius transforms O(n*2^n)
- Fast multi-point evaluation O(n*2^n)
- Subset convolution O(n^2*2^n)
- Williams SAT O*(2^{n/2})
- OV solvers: brute-force, packed, polynomial

## L6: Canonical Problems (Complete)

- Orthogonal Vectors (OV)
- CNF-SAT (Williams algorithm)
- MOD_m (modular counting)
- MAJORITY
- PARITY

## L7: Applications (Complete)

- Fine-grained complexity (OV -> conditional lower bounds)
- Circuit lower bounds (Razborov-Smolensky)
- Algorithm design (polynomial method for SAT)
- Boolean function analysis (influences, Fourier)
- Cryptography (algebraic attacks via ANF)
- Learning theory (low-degree algorithm)

## L8: Advanced Topics (Complete)

- Razborov-Smolensky method
- Williams ACC lower bound (NEXP not in ACC0)
- Natural proofs barrier (Razborov-Rudich)
- Algebraic circuit complexity
- Fourier concentration theory

## L9: Research Frontiers (Partial)

- Meta-complexity (MCSP) - implemented
- P vs NP barriers - documented
- GCT (Geometric Complexity Theory) - future
- Quantum polynomial method - future
