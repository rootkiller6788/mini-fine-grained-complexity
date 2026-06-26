# Knowledge Graph ? mini-orthogonal-vectors

## L1: Definitions (Complete)
- inary_vector_t ? Binary vector in {0,1}^d, packed word representation
- ector_set_t ? Collection of n vectors sharing dimension d
- ov_instance_t ? OV instance: two vector sets A and B
- ov_result_t ? Query result with indices and statistics
- density_stats_t ? Vector set density statistics
- Dot product <a,b> for binary vectors
- Orthogonality: <a,b> = 0
- Hamming weight, Hamming distance
- Jaccard similarity, Cosine similarity

## L2: Core Concepts (Complete)
- OV Conjecture: n^{2-o(1)} for d = ?(log n)
- Three regimes: d = o(log n), d = c?log n, d = ?(log n)
- SETH ? OV Conjecture (Williams 2005)
- Subquadratic equivalence class
- Fine-grained complexity framework

## L3: Mathematical Structures (Complete)
- Packed bit representation: ceil(d/32) words
- Word-level parallelism: AND, OR, XOR, NOT
- Population count (popcount)
- GF(2) vector space operations
- Gaussian elimination over GF(2)
- Bit frequency analysis

## L4: Fundamental Laws (Complete)
- Williams' Theorem: SETH ? OV Conjecture
- Williams' reduction: k-SAT ? OV
- OV lower bound for d = ?(log n)
- Sparsification lemma connection
- Polynomial method lower bound
- Communication complexity lower bound: ?(n)

## L5: Algorithms/Methods (Complete)
- Brute force: O(n??d)
- Brute force with early termination
- Sparse-first ordering heuristic
- Williams' algorithm: random projection, O(n^{2-1/O(c)})
- Meet-in-the-Middle: O(n?2^{d/2})
- Four-Russians speedup: O(n??d/log n)
- Matrix multiplication approach: O(n^?)
- Polynomial method (AWY 2015)
- Light-sparse optimization: O(n????d)
- Locality-Sensitive Hashing: bit-sampling
- Counting variants, batch operations
- Streaming algorithm

## L6: Canonical Problems (Complete)
- OV itself (the canonical problem)
- Pattern Matching with Wildcards (reduction)
- Graph Diameter in sparse graphs (Roditty-Williams 2013)
- Edit Distance / Levenshtein distance (Bringmann-Kunnemann 2015)
- Longest Common Subsequence (Abboud-Backurs-Williams 2015)
- Subset Sum (Abboud-Bringmann-Hermelin-Shabtay 2019)
- OV Equivalence Class documentation
- Frechet Distance, DTW (documented)

## L7: Applications (Complete ? 4 applications)
- Database anti-join optimization
- CRISPR off-target prediction (computational biology)
- Code similarity / plagiarism detection
- Collaborative filtering (diverse user detection)
- Parameter sweep / empirical analysis

## L8: Advanced Topics (Partial+ ? 5 topics)
- Dynamic OV data structure (insert, delete, query)
- Approximate OV (?-relaxed orthogonality)
- Polynomial method degree analysis
- Communication complexity simulation
- Instance difficulty estimation
- GF(2) rank computation
- Phase transition analysis

## L9: Research Frontiers (Partial ? documented)
- OV ? 3SUM conjecture implication (negative result documented)
- OV in non-binary domains (documented)
- Equivalence class ongoing research
- Implication chain visualization
