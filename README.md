# Mini Fine-Grained Complexity

A collection of **from-scratch, zero-dependency C implementations** of fine-grained complexity theory — the study of exact polynomial-time lower bounds for fundamental computational problems. Each module maps to MIT and Stanford advanced complexity courses, translating conditional lower bounds, equivalence classes, and hardness conjectures (SETH, ETH, OVC, 3SUM, APSP) into runnable C code.

## Sub-Modules

| Sub-Module | Topics | Key Courses |
|-----------|--------|-------------|
| [mini-3sum-apsp-conjectures](mini-3sum-apsp-conjectures/) | 3SUM conjecture, APSP conjecture, subcubic equivalence, min-plus product, reduction chains | MIT 6.8410, Stanford CS367 |
| [mini-conditional-lower-bounds](mini-conditional-lower-bounds/) | CLB framework, SETH/ETH/OVC/3SUM/APSP hypotheses, reduction web, string hardness | MIT 6.8410, Stanford CS254 |
| [mini-edit-distance-lcs-hardness](mini-edit-distance-lcs-hardness/) | Edit distance, LCS, Fréchet distance, sequence alignment, quadratic lower bounds under SETH | MIT 6.8410, Stanford CS367 |
| [mini-equivalence-classes](mini-equivalence-classes/) | Subcubic equivalence class, subquadratic equivalence class, 3SUM class, fine-grained reductions | MIT 6.8410, Stanford CS367 |
| [mini-k-clique-hardness](mini-k-clique-hardness/) | k-Clique parameterized hardness, ETH/SETH lower bounds, FPT algorithms, color-coding | MIT 6.8410, Stanford CS367 |
| [mini-orthogonal-vectors](mini-orthogonal-vectors/) | OV problem, SETH→OV reduction, OV→Edit/LCS/Diameter reductions, pattern matching | MIT 6.8410, Stanford CS367 |
| [mini-polynomial-method-consequences](mini-polynomial-method-consequences/) | Polynomial method, algebraic complexity, faster-than-brute-force via algebraic techniques | MIT 6.8420, Stanford CS254 |
| [mini-seth-strong-eth](mini-seth-strong-eth/) | SETH/ETH formalization, exponential complexity, CNF generation, sparsification | MIT 6.8410, Stanford CS254 |

## Design Philosophy

- **Zero external dependencies** — pure C (C99/C11), only `libc` and `libm`
- **Self-contained modules** — each directory has its own `include/`, `src/`, `tests/`, and `Makefile`
- **Theory-to-code mapping** — every module includes `docs/` with conjecture formalizations and course-alignment notes
- **Practical benchmarking** — empirical time-complexity measurement, scaling analysis, and reduction verification

## Building

Each module is standalone. Navigate to a module directory and run:

```bash
cd mini-3sum-apsp-conjectures
make all    # build everything
make test   # run tests
```

Requires **GCC** and **GNU Make**.

## Project Structure

```
mini-fine-grained-complexity/
├── mini-3sum-apsp-conjectures/         # 3SUM & APSP conjectures, subcubic equivalence
├── mini-conditional-lower-bounds/      # Conditional lower bounds framework, reduction web
├── mini-edit-distance-lcs-hardness/    # Edit distance, LCS, Fréchet distance hardness
├── mini-equivalence-classes/           # Subcubic, subquadratic, 3SUM equivalence classes
├── mini-k-clique-hardness/             # k-Clique parameterized complexity, ETH/SETH bounds
├── mini-orthogonal-vectors/            # OV problem, SETH connection, quadratic reductions
├── mini-polynomial-method-consequences/# Polynomial method in fine-grained complexity
└── mini-seth-strong-eth/               # SETH & ETH formalization, exponential complexity
```

## License

MIT
