# doc-1: SETH and ETH Fundamentals

## L1-L3: Core Definitions and Mathematical Structure

### 1.1 The Exponential Time Hypothesis (ETH)

**Statement (Impagliazzo-Paturi 1999):**
3-SAT cannot be solved in time 2^{o(n)} where n is the number of variables.

**Formal definition:**
There exists a constant ε₀ > 0 such that 3-SAT requires time Ω(2^{ε₀·n}) in the worst case on a deterministic Turing machine.

**Key implications:**
- P ≠ NP (since SAT ⊂ P would give polynomial time for 3-SAT)
- All NP-complete problems require exponential time
- The exact constant ε₀ for 3-SAT is unknown (best known: s₃ = 0.3863)

### 1.2 The Strong Exponential Time Hypothesis (SETH)

**Statement (Impagliazzo-Paturi-Zane 2001):**
For every ε > 0, there exists an integer k ≥ 3 such that k-SAT cannot be solved in time O(2^{(1-ε)n}) on instances with n variables.

**Equivalent formulation:**
lim_{k→∞} s_k = 1, where s_k = inf{c : k-SAT can be solved in O(2^{c·n}) time}

**Known s_k values (PPSZ 2005):**
- s₃ = 0.3863
- s₄ = 0.5548
- s₅ = 0.6500
- s₆ = 0.7162
- s₇ = 0.7627
- s₈ = 0.7977
- s₉ = 0.8245
- s₁₀ = 0.8453
- sₖ ≈ 1 - 1.5/k for large k

### 1.3 Relationship: ETH vs SETH

| Property | ETH | SETH |
|----------|-----|------|
| Formulated | 1999 | 2001 |
| Concerned with | 3-SAT specifically | All k-SAT |
| Lower bound | 2^{Ω(n)} | 2^{(1-o(1))n} |
| Implies | P ≠ NP | P ≠ NP |
| Refuted by | 2^{o(n)} 3-SAT algorithm | k-SAT algorithm with exponent < 1 for all k |

**Theorem (IPZ 2001):** ETH ⇒ SETH (via Sparsification Lemma)

### 1.4 CNF Formula Data Structures

**literal_t:** A variable index (1-based) with a negation bit
**clause_t:** An array of literals representing a disjunction
**cnf_formula_t:** An array of clauses representing a conjunction
**assignment_t:** An array of values {-1=unset, 0=false, 1=true}

### 1.5 Core Operations

- **Unit propagation:** If a clause has all but one literal falsified, force the remaining literal true
- **Pure literal elimination:** If a variable appears only positively or only negatively, assign it accordingly
- **Subsumption:** If clause A's literals are a subset of clause B's, then B is redundant

## References

1. Impagliazzo, R., Paturi, R. (1999). "Complexity of k-SAT." CCC 1999.
2. Impagliazzo, R., Paturi, R., Zane, F. (2001). "Which Problems Have Strongly Exponential Complexity?" JCSS 63(4):512-530.
3. Paturi, R., Pudlak, P., Saks, M.E., Zane, F. (2005). "An Improved Exponential-time Algorithm for k-SAT." JACM 52(3):337-364.
