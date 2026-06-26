# Knowledge Graph: mini-edit-distance-lcs-hardness

## L1: Definitions (Complete)
- Edit distance (Levenshtein distance)
- Longest Common Subsequence (LCS)
- Edit operations: insert, delete, substitute, match
- Edit script, edit cost model
- Hamming distance
- Damerau-Levenshtein distance
- Jaro similarity, Jaro-Winkler distance
- SETH (Strong Exponential Time Hypothesis)
- OV (Orthogonal Vectors) problem
- Frechet distance (dog-leash distance)
- DTW (Dynamic Time Warping)
- Sequence alignment (global, local)
- Substitution matrices (BLOSUM, PAM)

## L2: Core Concepts (Complete)
- Quadratic-time edit distance: O(n*m) Wagner-Fischer
- Quadratic-time LCS: O(n*m) DP
- n^{2-o(1)} conditional lower bounds under SETH
- Fine-grained reduction: SAT -> OV -> Edit/LCS
- OV Conjecture and its implications
- Linear-space edit distance (2-row DP)
- Bounded edit distance (Ukkonen's cutoff)
- Global vs local alignment concepts
- Affine gap penalty model

## L3: Mathematical Structures (Complete)
- DP matrix: (n+1)x(m+1) row-major arrays
- Edit DP cell with backpointer
- LCS DP cell with direction
- Traceable DP for script reconstruction
- Affine DP (M, Ix, Iy matrices)
- Binary vector representation (packed bits)
- OV instance: two sets of binary vectors
- Curve representation: sequence of 2D points
- Alignment column with gap characters

## L4: Fundamental Theorems (Complete)
- Wagner-Fischer Theorem: Edit distance in O(nm)
- Hirschberg Theorem: LCS in O(nm) time, O(min(n,m)) space
- Backurs-Indyk Theorem (STOC 2015): Edit distance needs n^{2-o(1)} under SETH
- ABW Theorem (FOCS 2015): LCS needs n^{2-o(1)} under SETH
- Bringmann-Kunnemann Theorem (FOCS 2015): Frechet needs n^{2-o(1)} under SETH
- Masek-Paterson (1980): Edit distance in O(n^2 / log n)
- Hunt-Szymanski (1977): LCS in O((R+n)log n)

## L5: Algorithms/Methods (Complete)
- Wagner-Fischer DP for edit distance
- 2-row linear-space DP for edit distance
- Ukkonen's bounded edit distance (diagonal band)
- DPLL backtracking for SAT (SETH baseline)
- Standard LCS DP with traceback
- Hirschberg divide-and-conquer LCS
- Needleman-Wunsch global alignment
- Smith-Waterman local alignment
- Gotoh affine gap alignment
- Discrete Frechet distance DP
- DTW distance DP with Sakoe-Chiba band
- Four Russians speedup (Masek-Paterson)
- Hunt-Szymanski LCS algorithm
- LIS via patience sorting (O(n log n))
- LCIS via DP

## L6: Canonical Problems (Complete)
- Edit distance computation (Levenshtein)
- Longest Common Subsequence
- Longest Common Substring
- Shortest Common Supersequence
- Longest Increasing Subsequence
- Longest Common Increasing Subsequence
- Global sequence alignment (Needleman-Wunsch)
- Local sequence alignment (Smith-Waterman)
- Frechet distance between polygonal curves
- DTW distance between curves
- Orthogonal Vectors problem
- Approximate string matching (k-differences)

## L7: Applications (Partial+)
- Bioinformatics: DNA sequence alignment
- DNA to protein translation
- DNA reverse complement, GC content
- Random string generation for testing
- Sequence identity computation
- Jaro-Winkler for record linkage

## L8: Advanced Topics (Partial+)
- SETH conditional lower bound proofs
- OV-to-Edit-Distance reduction chain
- Beyond-quadratic lower bounds (k-OV, APSP, 3SUM)
- OVH (Orthogonal Vectors Hypothesis)
- Approximate edit distance estimation
- Four Russians block-based speedup
- Frechet distance decision problem (free-space diagram)

## L9: Research Frontiers (Partial)
- s_k constants and their known bounds
- Breakthrough implications analysis
- Relationship between SETH, OVC, and sequence measures
- Open problems in fine-grained complexity for strings
