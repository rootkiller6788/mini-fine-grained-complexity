# Gap Report: mini-edit-distance-lcs-hardness

## Missing Items (by priority)

### Priority 1 (Critical for Complete)
- None. All critical items are implemented.

### Priority 2 (Should Have)
- Full BLOSUM62 matrix data (currently simplified to identity matrix)
- True Hirschberg linear-space edit script reconstruction (currently uses full DP internally)
- Affine gap alignment with proper 3-matrix DP (currently delegates to linear-gap NW)
- Four Russians with boundary state propagation (currently simplified block DP)

### Priority 3 (Nice to Have)
- Sublinear-space edit distance approximation (Andoni-Krauthgamer-Onak full implementation)
- True Hunt-Szymanski algorithm for sparse LCS (currently delegates to standard LCS)
- Continuous Frechet distance exact algorithm (currently discrete approximation)
- More substitution matrices (BLOSUM62, PAM250, etc.)
- Progressive MSA with guide tree (UPGMA/neighbor-joining)

## L9 Research Frontiers Gaps
- No quantum algorithms for edit distance
- No streaming algorithms for edit distance
- No implementation of recent breakthroughs (e.g., O(n^{2-2/7}) for binary alphabet)
- No implementation of metric embeddings for edit distance

## Summary
The module provides comprehensive coverage of the edit distance and LCS hardness landscape. All core theorems, algorithms, and canonical problems are implemented and tested. Advanced topics and research frontiers are documented with partial implementations where feasible.
