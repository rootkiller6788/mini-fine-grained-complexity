# doc-2: Sparsification Lemma and Isolation Lemma

## L4: Fundamental Theorems

### 2.1 Sparsification Lemma (IPZ 2001)

**Statement:**
For any k ≥ 3 and ε > 0, any k-CNF formula F with n variables can be expressed as:

    F = F₁ ∨ F₂ ∨ ... ∨ F_t

where t ≤ 2^{εn} and each F_i is a k-CNF formula with at most C(k, ε) · n clauses, where C(k, ε) is a constant depending only on k and ε.

**Proof Technique:**
1. **Sunflower Lemma (Erdős-Rado 1960):** A sunflower is a collection of sets S₁, ..., S_p where S_i ∩ S_j = C (the "core") for all i ≠ j, and the "petals" P_i = S_i \ C are pairwise disjoint.
2. If a family of k-clauses contains a sunflower with p petals, we can replace it with just the core clause C, reducing the number of clauses.
3. **Encoding:** View each clause as a partial assignment. The sunflower lemma guarantees that a large enough clause family contains a sunflower.
4. **Branching:** Branch on subsets of variables to identify the sunflowers.

**C(k, ε) bound:** C(k, ε) = (k/ε)^{O(k)} — doubly exponential, not practical for SAT solving but sufficient for theory.

**Key application: ETH ⇒ SETH**
- If SETH is false (∃δ > 0: ∀k, k-SAT ∈ TIME(2^{(1-δ)n})), then sparsification with ε = δ/2 followed by k-SAT solving gives a 2^{(1-δ/2)n} algorithm for 3-SAT, contradicting ETH.

### 2.2 Isolation Lemma (Valiant-Vazirani 1986)

**Statement:**
For any CNF formula F on n variables, randomly hashing with O(log n) bits produces a formula that has exactly one satisfying assignment (if F is satisfiable) with probability ≥ 1/(4n).

**Application:**
- Derandomization of the PPSZ algorithm
- Reducing SAT to Unique-SAT
- Randomized polynomial hierarchy collapse: NP ⊆ RP^{PromiseUP}

**Implementation:**
```c
isolation_result_t *isolation_apply(const cnf_formula_t *f, uint64_t seed);
```
- Adds random constraints using O(log n) hash bits
- Success probability ≥ 1/(4n)
- Seed determines the random hash function

### 2.3 Duality Between Clause Width and Time

**Calabro-Impagliazzo-Paturi (2006):**
There is a fundamental duality between the clause width k and the running time exponent s_k. The Sparsification Lemma formalizes this: reducing clause width (making formulas sparser) reduces the time needed to solve k-SAT.

**Key insight:**
- Dense formulas (high clause-to-variable ratio): Harder to find satisfying assignment, but easier to detect unsatisfiability
- Sparse formulas: Easier to find satisfying assignment via local search
- Phase transition at α_c ≈ 2^k · ln(2) separates these regimes

## References

1. Impagliazzo, R., Paturi, R., Zane, F. (2001). JCSS 63(4):512-530.
2. Valiant, L.G., Vazirani, V.V. (1986). "NP is as easy as detecting unique solutions." TCS 47:85-93.
3. Calabro, C., Impagliazzo, R., Paturi, R. (2006). "A Duality between Clause Width and Time for SAT." CCC 2006.
