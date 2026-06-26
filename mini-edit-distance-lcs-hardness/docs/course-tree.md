# Course Tree: mini-edit-distance-lcs-hardness

## Prerequisites
- mini-seth-strong-eth: SETH definitions, k-SAT, exponential complexity
- mini-orthogonal-vectors: OV problem, reductions from SAT
- mini-k-clique-hardness: Fine-grained reduction techniques
- Basic algorithms: Dynamic programming, string algorithms
- Basic complexity: P vs NP, polynomial-time reductions

## Dependencies (what depends on this module)
- mini-conditional-lower-bounds: Edit/LCS as case studies
- mini-equivalence-classes: Edit/LCS in fine-grained equivalence classes
- mini-polynomial-method-consequences: Implications of quadratic lower bounds

## Knowledge Flow
```
SETH (mini-seth-strong-eth)
  |
  v
OV Conjecture (mini-orthogonal-vectors)
  |
  +---> Edit Distance Hardness [THIS MODULE]
  |       |
  |       +---> Backurs-Indyk (2015): SAT -> OV -> Edit Distance
  |       +---> ABW (2015): SAT -> OV -> LCS
  |       +---> Bringmann-Kunnemann (2015): SAT -> OV -> Frechet
  |
  +---> Conditional Lower Bounds (mini-conditional-lower-bounds)
  +---> Equivalence Classes (mini-equivalence-classes)
  +---> Polynomial Method (mini-polynomial-method-consequences)
```

## Key Papers Implemented
1. Levenshtein (1966) - Binary codes for deletions/insertions
2. Wagner & Fischer (1974) - String-to-string correction
3. Hirschberg (1975) - Linear space LCS
4. Needleman & Wunsch (1970) - Global alignment
5. Smith & Waterman (1981) - Local alignment
6. Masek & Paterson (1980) - Four Russians for edit distance
7. Ukkonen (1985) - Approximate string matching
8. Backurs & Indyk (STOC 2015) - Edit distance hardness under SETH
9. Abboud, Backurs, Williams (FOCS 2015) - LCS hardness under SETH
10. Bringmann & Kunnemann (FOCS 2015) - Frechet/DTW hardness under SETH
