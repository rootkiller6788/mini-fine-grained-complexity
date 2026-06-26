# mini-edit-distance-lcs-hardness

## Status: COMPLETE ✅

**Edit Distance and LCS quadratic lower bounds under SETH.**

- Backurs-Indyk (2015): Edit Distance needs n^{2-o(1)} under SETH
- ABW (2015): LCS needs n^{2-o(1)} under SETH
- Bringmann-Kunnemann (2015): Frechet distance needs n^{2-o(1)} under SETH

## Knowledge Coverage

| Level | Name | Status |
|-------|------|--------|
| L1 | Definitions | Complete |
| L2 | Core Concepts | Complete |
| L3 | Mathematical Structures | Complete |
| L4 | Fundamental Theorems | Complete |
| L5 | Algorithms/Methods | Complete |
| L6 | Canonical Problems | Complete |
| L7 | Applications | Complete |
| L8 | Advanced Topics | Complete |
| L9 | Research Frontiers | Partial |

## Core Definitions (L1)
- `edit_op_t`, `edit_script_t`, `edit_cost_model_t`: Edit operations and costs
- `lcs_result_t`, `lcs_match_t`: LCS result with matching positions
- `alignment_result_t`, `alignment_scoring_t`: Sequence alignment
- `curve_t`, `point2d_t`: Curves for Frechet distance
- `binary_vector_t`, `ov_instance_t`: Orthogonal Vectors problem

## Core Theorems (L4)
- **Wagner-Fischer (1974)**: Edit distance in O(nm) DP
- **Hirschberg (1975)**: LCS in O(nm) time, O(min(n,m)) space
- **Backurs-Indyk (STOC 2015)**: ED requires n^{2-o(1)} under SETH
- **ABW (FOCS 2015)**: LCS requires n^{2-o(1)} under SETH
- **Bringmann-Kunnemann (FOCS 2015)**: Frechet/DTW requires n^{2-o(1)} under SETH
- **Masek-Paterson (1980)**: ED in O(n^2/log n) via Four Russians

## Core Algorithms (L5)
- Wagner-Fischer DP: O(nm) edit distance
- 2-row linear-space DP: O(min(n,m)) space
- Ukkonen bounded DP: O(k*min(n,m))
- Hirschberg divide-and-conquer LCS
- Needleman-Wunsch global alignment
- Smith-Waterman local alignment
- Discrete Frechet distance DP
- DTW with Sakoe-Chiba band
- Four Russians block-based speedup
- LIS via patience sorting

## Canonical Problems (L6)
- Edit distance (Levenshtein, Damerau-Levenshtein)
- Longest Common Subsequence / Substring
- Shortest Common Supersequence
- Global/Local sequence alignment
- Frechet distance, DTW
- Orthogonal Vectors (OV)

## Course Mapping
- MIT 6.890: Fine-Grained Complexity
- Stanford CS254: Computational Complexity
- Berkeley CS278: Computational Complexity
- CMU 15-751: Theorist's Toolkit
- Princeton COS 522: Computational Complexity

## Build & Test
```
make          # build library (libedlcs.a)
make test     # run 55 tests (all passing)
make main     # build demo binary
make run      # run demo
make clean    # clean build artifacts
make examples # build example programs
```

## File Structure
```
include/         6 header files (842 lines)
  edit_distance.h, lcs.h, seth_hardness.h, alignment.h,
  frechet_distance.h, string_utils.h
src/             7 source files (2724 lines)
  edit_distance.c, lcs.c, seth_hardness.c, alignment.c,
  frechet_distance.c, string_utils.c, main.c
tests/           1 test file (193 lines, 55 tests)
examples/        3 example programs
docs/            5 knowledge documents
```
