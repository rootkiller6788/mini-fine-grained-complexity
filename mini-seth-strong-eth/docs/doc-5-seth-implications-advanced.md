# doc-5: SETH Implications and Advanced Topics

## L8-L9: Circuit Lower Bounds, Parameterized Complexity, and Empirical Validation

### 5.1 SETH and Circuit Complexity

**Monotone Circuit Lower Bounds:**
Under SETH, the monotone circuit complexity of k-SAT grows as n^{k/(k-1)}. This connects fine-grained time lower bounds to circuit size lower bounds.

**Function:** `seth_monotone_circuit_bound(n, k)` → n^{k/(k-1)}

**NC vs P Separation:**
SETH implies that SAT is not in NC: polynomial-size, polylog-depth circuits cannot solve SAT. The relationship is parameterized by k (clause width) and depth:
- `seth_implies_nc_sat_separation(k, depth)` → depth < k/2 ⇒ separation

### 5.2 ETH and Parameterized Complexity

**The W-Hierarchy and ETH:**

ETH provides **quantitative** lower bounds for FPT (Fixed-Parameter Tractable) problems, going beyond what W-hardness alone provides.

| Problem | W-class | ETH Lower Bound | Best Algorithm |
|---------|---------|-----------------|----------------|
| k-Clique | W[1]-complete | n^{Ω(k)} | O(n^{ωk/3}) where ω < 2.373 |
| k-Dominating Set | W[2]-complete | n^{Ω(k)} | O(n^{k+1}) |
| Hitting Set(k) | W[2]-complete | n^{Ω(k log k)} | n^{O(k)} |

**Key consequence:** If ETH holds, then FPT ≠ W[1]. This is stronger than P ≠ NP, showing that even with arbitrary parameter-dependence f(k), k-Clique needs super-polynomial time in n.

### 5.3 Treewidth and SETH

**Courcelle's Theorem** gives FPT algorithms parameterized by treewidth tw: O(f(tw) · n). Under SETH, the function f(tw) must be at least exponential:

| Problem | SETH-optimal base | Best known |
|---------|-------------------|------------|
| Vertex Cover | 2^{tw} | 2^{tw} · n |
| Dominating Set | 3^{tw} | 3^{tw} · n |
| Independent Set | 2^{tw} | 2^{tw} · n |
| 3-Coloring | 3^{tw} | 3^{tw} · n |
| Hamiltonicity | 4^{tw} | 4^{tw} · n |

Achieving a smaller base for any of these would refute SETH.

### 5.4 Set Partition and Subset Sum

**Set Partition** under SETH: requires 2^{(1-o(1))n} time.
- Best known: O(2^{n/2}) via meet-in-the-middle (Horowitz-Sahni 1974)
- SETH says this is essentially optimal

**Subset Sum** under SETH: O(2^{n/2}) is optimal.
- The meet-in-the-middle algorithm is believed to be optimal under SETH

### 5.5 NP-Intermediate Problems under ETH

ETH provides a finer lens on problems between P and NP-complete:

**Graph Isomorphism (GI):**
- Status: Not known NP-complete under ETH (unless PH collapses)
- Babai (2016): GI ∈ DTIME(2^{O((log n)^3)}) — quasipolynomial time
- ETH predicts: GI is not NP-complete but also not in P

**Integer Factoring:**
- Status: BQP (Shor's algorithm) but not known in P
- ETH predicts: Factoring is not NP-complete (otherwise cryptography collapses)

### 5.6 Empirical Validation of SETH

**Methodology:**
1. Generate random k-SAT instances at the phase transition
2. Run SAT solvers and measure effective exponent c in T(n) ≈ 2^{c·n}
3. Observe c as k increases
4. SETH prediction: c → 1 as k → ∞

**Implementation:**
```c
void check_seth_consistency(solver, n, k_min, k_max, num_instances);
seth_evidence_level_t assess_seth_evidence(measured_s_k, k_min, k_max);
extrapolation_result_t extrapolate_runtime(solver, k, n_min, n_max, n_step);
```

**Evidence Levels:**
- SETH_SUPPORTED: measured s_k increases monotonically to ~1
- SETH_WEAK_EVIDENCE: insufficient data
- SETH_TENSION: s_k plateaus below 0.9
- SETH_REFUTED: algorithm achieves exponent < s_k - 0.1

### 5.7 k-SAT Phase Transition and SETH

The phase transition at clause-to-variable ratio α_c ≈ 2^k · ln 2 is where:
- Below α_c: formulas almost surely satisfiable (easy region)
- Above α_c: formulas almost surely unsatisfiable (easy region)
- Near α_c: exponentially hard instances

SETH implies that near-threshold instances truly require exponential time for large k.

### 5.8 Consequences if SETH is False

| Consequence | Impact |
|-------------|--------|
| P ≠ NP still possible | ETH may still hold |
| Faster O(N^{1.9}) Edit Distance | Huge practical impact for computational biology |
| Faster O(N^{1.9}) LCS | Impact for bioinformatics, diff tools |
| Sub-quadratic 3SUM | Impact for computational geometry |
| k-Dom Set in n^{0.9k} | Parameterized complexity revolution |

### 5.9 Open Problems

1. **Prove or refute SETH.** The central open problem. A single algorithm achieving exponent 1-ε for all k would refute SETH.
2. **Tighten s₃.** Is s₃ = 0.3863 optimal? Could there be a 2^{0.2n} algorithm for 3-SAT?
3. **Deterministic PPSZ.** Can PPSZ be derandomized without loss in the exponent?
4. **ETH without SETH.** Is it possible that ETH holds but SETH is false? Requires ETH ⇒ SETH proof to be nontrivial.
5. **Sparsification constants.** Improve C(k, ε) from (k/ε)^{O(k)} to something practical.

## References

1. Impagliazzo, R., Paturi, R., Zane, F. (2001). JCSS 63(4):512-530.
2. Williams, R. (2014). "Faster all-pairs shortest paths via circuit complexity." STOC 2014.
3. Lokshtanov, D., Marx, D., Saurabh, S. (2011). "Lower bounds based on ETH." Bulletin of the EATCS 105:41-72.
4. Cygan, M., et al. (2015). "Parameterized Algorithms." Springer.
5. Babai, L. (2016). "Graph Isomorphism in Quasipolynomial Time." arXiv:1512.03547.
