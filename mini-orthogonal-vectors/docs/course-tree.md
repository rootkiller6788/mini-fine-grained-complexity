# Course Tree ? mini-orthogonal-vectors

## Prerequisites (Knowledge Dependencies)

`
mini-complexity-foundations
??? mini-cook-levin-theorem          [L4: SAT is NP-complete]
??? mini-p-np-np-completeness        [L1: P, NP definitions]
??? mini-time-hierarchy-theorem      [L4: Time hierarchy]
??? mini-reductions-completeness     [L1: Reduction concept]

mini-fine-grained-complexity
??? mini-seth-strong-eth             [L2: SETH definition]
?   ??? mini-orthogonal-vectors      [THIS MODULE]
?       ??? L1: OV definitions
?       ??? L2: OV conjecture ? SETH
?       ??? L4: Williams' theorem
?       ??? L5: Williams' algorithm
?       ??? L6: Equivalence class
?       ??? L7-L9: Applications & frontiers
??? mini-k-clique-hardness           [depends on OV equivalence]
??? mini-3sum-apsp-conjectures       [parallel conjecture]
??? mini-edit-distance-lcs-hardness  [depends on OV ? ED reduction]
??? mini-conditional-lower-bounds    [methodology, depends on OV]
??? mini-equivalence-classes         [depends on OV class]
??? mini-polynomial-method-consequences [depends on polynomial method]
`

## Knowledge Flow Within This Module

`
L1: Definitions (binary_vector_t, ov_instance_t)
  ?
L3: Mathematical Structures (packed bits, GF(2) algebra)
  ?
L2: Core Concepts (OV conjecture, SETH connection)
  ?
L4: Fundamental Laws (Williams' theorem, lower bounds)
  ?
L5: Algorithms (brute force ? Williams ? 4-Russians ? MM ? LSH)
  ?
L6: Canonical Problems (reductions to PM, ED, LCS, Diameter, SubsetSum)
  ?
L7: Applications (database, CRISPR, code similarity, recommendation)
  ?
L8: Advanced Topics (dynamic OV, approximate OV, communication complexity)
  ?
L9: Research Frontiers (3SUM implication, non-binary OV)
`

## Downstream Dependencies

`
mini-orthogonal-vectors (THIS)
  ?
mini-edit-distance-lcs-hardness    [uses OV ? ED/LCS reductions]
mini-equivalence-classes           [uses OV equivalence class]
mini-conditional-lower-bounds      [uses OV as canonical hard problem]
mini-polynomial-method-consequences [uses polynomial method technique]
`

## Key Papers in Dependency Order
1. Williams (2005) ? SETH ? OV, Williams' algorithm
2. Roditty-Williams (2013) ? OV ? Diameter
3. Abboud-Williams-Yu (2015) ? Polynomial method, equivalence class
4. Bringmann-Kunnemann (2015) ? OV ? Edit Distance
5. Abboud-Backurs-Williams (2015) ? OV ? LCS, DTW
6. Abboud-Bringmann-Hermelin-Shabtay (2019) ? OV ? Subset Sum
