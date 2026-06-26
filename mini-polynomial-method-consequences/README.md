# mini-polynomial-method-consequences

**Polynomial Method Applications to Fine-Grained Complexity**

## Module Status: COMPLETE ✅

- L1-L6: Complete
- L7: Complete (6 applications)
- L8: Complete (Razborov-Smolensky, ACC lower bounds, natural proofs)
- L9: Partial (documented, MCSP implemented)

## Core Definitions (L1)

- `polynomial_t`: Multivariate polynomial over GF(2) or reals (sparse representation)
- `monomial_t`: Product of variables as bitmask
- `term_t`: Coefficient × monomial
- `gf2_t` / `gf_p_t`: Finite field elements and arithmetic
- `algebraic_circuit_t`: DAG computing polynomials (+ and × gates)
- `ov_instance_t`: Orthogonal Vectors problem instance
- `williams_cnf_t`: CNF formula for Williams' SAT algorithm

## Core Theorems (L4)

- **Zhegalkin's Theorem (1927)**: Every Boolean function f:{0,1}^n → {0,1} has a UNIQUE multi-linear polynomial over GF(2)
- **Fourier Uniqueness**: Unique spectral representation f(x) = Σ f̂(S)·χ_S(x)
- **Razborov's Theorem (1987)**: MAJORITY not in AC^0[p]; every AC^0 circuit of size s and depth d approximated by degree (log s)^{O(d)} polynomial over GF(p)
- **Smolensky's Theorem (1987)**: MOD_q not in AC^0[p] for distinct primes p,q; MOD_3 requires degree Ω(n) over GF(2)
- **Williams' Algorithm (2010/2014)**: CNF-SAT in O(2^{n-n/O(log m)}) via polynomial method; NEXP not in ACC^0

## Core Algorithms (L5)

- **Fast Walsh-Hadamard Transform**: O(n·2^n) Fourier transform on Boolean cube
- **Zeta/Möbius Transforms**: O(n·2^n) subset sum and inversion
- **Fast Multi-Point Boolean Evaluation**: Evaluate polynomial at ALL 2^n points in O(n·2^n)
- **Subset Convolution**: (f·g) via ranked Möbius, O(n²·2^n)
- **Williams SAT Solver**: O*(2^{n/2}) via split-and-polynomial
- **OV Solvers**: Brute-force O(n²·d), Packed O(n²), Polynomial O(n²·2^d)

## Key Numbers (L6)

| Function | GF(2)-degree | GF(3)-degree | Significance |
|----------|-------------|-------------|-------------|
| PARITY (MOD_2) | 1 | n | MOD_2 easy over GF(2), hard over GF(3) |
| MOD_3 | Ω(n) | ≤ 2 | MOD_3 hard over GF(2) ⇒ not in AC^0[2] |
| MAJORITY | Ω(n) | Ω(n) | Not in AC^0[p] for any p (Razborov) |
| AND_n | n | n | Degree equals n for exact representation |

## Course Mapping

- MIT 6.841/6.890: Circuit Complexity, Fine-Grained Complexity
- Stanford CS254/CS358: Lower Bounds, Circuit Complexity
- Berkeley CS278: Computational Complexity
- CMU 15-855: Graduate Complexity Theory
- Princeton COS 522/551: Computational/Advanced Complexity
- Cambridge Part III: Advanced Complexity
- Oxford/ETH: Advanced Complexity Theory

## Build & Test

```
make        # build library, tests, examples, benches
make test   # run 40 tests (all passing)
make run-examples  # run 3 examples
make run-benches   # run 2 benchmarks
make clean  # clean build artifacts
```

## Module Structure

```
mini-polynomial-method-consequences/
├── Makefile
├── README.md
├── include/poly_method.h      # All type definitions and declarations
├── src/
│   ├── poly_method.c          # Core polynomial operations & transforms
│   ├── ov_poly.c              # OV problem via polynomial method
│   ├── williams.c             # Williams' SAT algorithm
│   ├── circuit_poly.c         # Circuit complexity connections
│   ├── lower_bounds.c         # Lower bounds via polynomial representations
│   └── main.c                 # Demo entry point
├── tests/test_main.c          # 40 comprehensive tests
├── examples/
│   ├── ov_example.c           # OV polynomial method demonstration
│   ├── williams_example.c     # Williams SAT algorithm demonstration
│   └── circuit_example.c      # Circuit lower bounds demonstration
├── benches/
│   ├── ov_bench.c             # OV solver performance benchmarks
│   └── williams_bench.c       # Williams SAT algorithm benchmarks
└── docs/
    ├── knowledge-graph.md     # L1-L9 knowledge coverage
    ├── coverage-report.md     # Per-level completion assessment
    ├── gap-report.md          # Missing items and priorities
    ├── course-alignment.md    # 9-school curriculum mapping
    └── course-tree.md         # Prerequisite dependencies
```
