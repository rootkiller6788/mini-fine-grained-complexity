#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <math.h>

/* ============================================================================
 * mini-edit-distance-lcs-hardness: Sequence Alignment
 *
 * Global and local sequence alignment algorithms, including
 * Needleman-Wunsch, Smith-Waterman, and affine gap penalty models.
 *
 * Reference: Needleman & Wunsch (1970) JMB 48(3):443-453
 *            Smith & Waterman (1981) JMB 147(1):195-197
 *            Gotoh (1982) JMB 162(3):705-708 (affine gaps)
 *
 * L1: Alignment definitions, scoring schemes
 * L2: Global vs local alignment concepts
 * L3: Substitution matrices (BLOSUM, PAM-like)
 * L5: Needleman-Wunsch, Smith-Waterman algorithms
 * L7: Bioinformatics applications
 * ============================================================================ */

/* ---- L1: Alignment Type Definitions ---- */

/* Alignment types */
typedef enum {
    ALIGN_GLOBAL = 0,     /* Needleman-Wunsch: entire sequences */
    ALIGN_LOCAL  = 1,     /* Smith-Waterman: best subsequence match */
    ALIGN_SEMIGLOBAL = 2, /* End-gap free: no penalty for terminal gaps */
    ALIGN_OVERLAP = 3     /* Overlap alignment */
} alignment_type_t;

/* A single aligned column in the alignment output */
typedef struct {
    char    char_a;     /* Character from sequence A (or '-' for gap) */
    char    char_b;     /* Character from sequence B (or '-' for gap) */
    int32_t score;      /* Score contribution of this column */
} alignment_column_t;

/* Alignment result */
typedef struct {
    alignment_column_t *columns;
    int32_t             num_columns;
    int32_t             total_score;
    int32_t             identities;    /* Number of exact matches */
    int32_t             gaps;          /* Total number of gap positions */
    double              percent_identity;
    char               *aligned_a;     /* Aligned string A with gaps */
    char               *aligned_b;     /* Aligned string B with gaps */
} alignment_result_t;

/* ---- L2: Scoring Models ---- */

/* Simple scoring scheme: match/mismatch/gap */
typedef struct {
    int32_t match_score;
    int32_t mismatch_penalty;
    int32_t gap_open_penalty;
    int32_t gap_extend_penalty;
    bool    use_affine_gaps;  /* If true, use affine gap model */
} alignment_scoring_t;

/* A substitution matrix entry */
typedef struct {
    char    from;
    char    to;
    int32_t score;
} sub_matrix_entry_t;

/* Substitution matrix (e.g., BLOSUM62, PAM250) */
typedef struct {
    sub_matrix_entry_t *entries;
    int32_t             num_entries;
    char               *alphabet;
    int32_t             alphabet_size;
    int32_t            *score_matrix;  /* Flattened alphabet_size x alphabet_size */
    int32_t             default_score; /* For characters not in alphabet */
} substitution_matrix_t;

/* ---- L2: Alignment Algorithms ---- */

/* Needleman-Wunsch global alignment with linear gap penalty.
 * O(n*m) time, O(n*m) space.
 * Returns the optimal global alignment. */
alignment_result_t *align_needleman_wunsch(const char *a, const char *b,
                                            const alignment_scoring_t *scoring);

/* Smith-Waterman local alignment.
 * O(n*m) time, O(n*m) space.
 * Finds the best-scoring local alignment (subsequence match). */
alignment_result_t *align_smith_waterman(const char *a, const char *b,
                                          const alignment_scoring_t *scoring);

/* Needleman-Wunsch with affine gap penalties (Gotoh 1982).
 * Uses three DP matrices: M (match), I_x (insert), I_y (delete).
 * O(n*m) time, O(n*m) space. */
alignment_result_t *align_affine_gap(const char *a, const char *b,
                                      const alignment_scoring_t *scoring);

/* Semiglobal alignment: free end gaps.
 * Useful for aligning a short sequence to a long one. */
alignment_result_t *align_semiglobal(const char *short_seq,
                                      const char *long_seq,
                                      const alignment_scoring_t *scoring);

/* Overlap alignment: find best overlap between suffix of A and prefix of B */
alignment_result_t *align_overlap(const char *a, const char *b,
                                   const alignment_scoring_t *scoring);

/* ---- Substitution Matrices ---- */

/* Load a built-in substitution matrix.
 * Available: "blosum62", "blosum50", "pam250", "identity", "dna_simple"
 * Returns NULL if name not recognized. */
substitution_matrix_t *substitution_matrix_load(const char *name);
void substitution_matrix_destroy(substitution_matrix_t *mat);
int32_t substitution_score(const substitution_matrix_t *mat, char a, char b);

/* ---- Alignment Utilities ---- */

/* Free an alignment result */
void alignment_result_destroy(alignment_result_t *r);

/* Print human-readable alignment */
void alignment_result_print(const alignment_result_t *r);

/* Compute alignment statistics */
double alignment_percent_identity(const alignment_result_t *r);
double alignment_percent_gaps(const alignment_result_t *r);
double alignment_e_value(int32_t score, int32_t m, int32_t n,
                          double lambda, double k);

/* ---- L3: DP Table for Alignment ---- */

typedef struct {
    int32_t *M;    /* Match/mismatch matrix */
    int32_t *Ix;   /* Insert in X (gap in Y) */
    int32_t *Iy;   /* Insert in Y (gap in X) */
    int32_t  rows;
    int32_t  cols;
} affine_dp_table_t;

affine_dp_table_t *affine_dp_create(int32_t len_a, int32_t len_b);
void affine_dp_destroy(affine_dp_table_t *dp);
void affine_dp_fill(affine_dp_table_t *dp, const char *a, const char *b,
                     const alignment_scoring_t *scoring);

/* ---- L5: Multiple Sequence Alignment (MSA) ---- */

/* Progressive multiple alignment using guide tree.
 * Uses pairwise alignment to build up a multiple alignment.
 * O(k^2 * n^2) for k sequences of length ~n. */
alignment_result_t *align_multiple_progressive(char **sequences,
                                                int32_t num_seqs);

/* Consensus sequence from multiple alignment */
char *alignment_consensus(const alignment_result_t *msa);

/* ---- L7: Bioinformatics Utilities ---- */

/* Translate DNA to protein (standard genetic code).
 * Input DNA sequence, output amino acid sequence.
 * Both must be null-terminated. */
char *dna_to_protein(const char *dna);

/* Reverse complement of a DNA sequence */
char *dna_reverse_complement(const char *dna);

/* GC content of a DNA sequence (0.0 to 1.0) */
double dna_gc_content(const char *dna);

/* Compute optimal alignment score for scoring-based alignments */
int32_t align_score_only(const char *a, const char *b,
                          const alignment_scoring_t *scoring,
                          alignment_type_t type);

/* Pairwise sequence identity: number of exact matches / aligned length */
double pairwise_identity(const char *a, const char *b);

#endif /* ALIGNMENT_H */
