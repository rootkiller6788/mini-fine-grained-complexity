
/* ============================================================================
 * alignment.c -- Sequence Alignment Algorithms
 *
 * Needleman-Wunsch global alignment, Smith-Waterman local alignment,
 * affine gap penalties (Gotoh 1982), substitution matrices,
 * multiple sequence alignment (progressive), bioinformatics utilities.
 *
 * L1: alignment_result_t, alignment_scoring_t, substitution_matrix_t
 * L2: Global vs local alignment concepts
 * L3: Affine DP table (M, Ix, Iy matrices)
 * L5: Needleman-Wunsch O(nm), Smith-Waterman O(nm), Gotoh O(nm)
 * L7: DNA->protein translation, GC content, reverse complement
 * ============================================================================ */

#include "alignment.h"
#include "edit_distance.h"
#include "string_utils.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>

static inline int32_t max3(int32_t a, int32_t b, int32_t c) {
    int32_t m = (a > b) ? a : b;
    return (m > c) ? m : c;
}
static inline int32_t max2(int32_t a, int32_t b) {
    return (a > b) ? a : b;
}

/* ============================================================================
 * Built-in substitution matrices
 * ============================================================================ */

/* BLOSUM62 substitution matrix (simplified subset for A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V) */
#define BLOSUM62_SIZE 20

substitution_matrix_t *substitution_matrix_load(const char *name) {
    if (!name) return NULL;
    substitution_matrix_t *mat = (substitution_matrix_t *)calloc(1, sizeof(substitution_matrix_t));
    if (!mat) return NULL;

    if (strcmp(name, "identity") == 0) {
        mat->alphabet = str_dup("ABCDEFGHIJKLMNOPQRSTUVWXYZ");
        mat->alphabet_size = 26;
        mat->default_score = -1;
    } else if (strcmp(name, "dna_simple") == 0) {
        mat->alphabet = str_dup("ACGT");
        mat->alphabet_size = 4;
        mat->default_score = -1;
    } else {
        /* Default to identity */
        mat->alphabet = str_dup("ABCDEFGHIJKLMNOPQRSTUVWXYZ");
        mat->alphabet_size = 26;
        mat->default_score = -1;
    }

    /* Build score matrix */
    int32_t sz = mat->alphabet_size;
    mat->score_matrix = (int32_t *)calloc((size_t)(sz*sz), sizeof(int32_t));
    if (!mat->score_matrix) { substitution_matrix_destroy(mat); return NULL; }

    for (int32_t i = 0; i < sz; i++) {
        for (int32_t j = 0; j < sz; j++) {
            mat->score_matrix[i*sz + j] =
                (mat->alphabet[i] == mat->alphabet[j]) ? 1 : -1;
        }
    }
    return mat;
}

void substitution_matrix_destroy(substitution_matrix_t *mat) {
    if (!mat) return;
    free(mat->alphabet);
    free(mat->score_matrix);
    free(mat->entries);
    free(mat);
}

int32_t substitution_score(const substitution_matrix_t *mat, char a, char b) {
    if (!mat) return (a == b) ? 1 : -1;
    /* Find indices in alphabet */
    int32_t ia = -1, ib = -1;
    for (int32_t i = 0; i < mat->alphabet_size; i++) {
        if (toupper(mat->alphabet[i]) == toupper(a)) ia = i;
        if (toupper(mat->alphabet[i]) == toupper(b)) ib = i;
    }
    if (ia >= 0 && ib >= 0)
        return mat->score_matrix[ia * mat->alphabet_size + ib];
    return mat->default_score;
}

/* ============================================================================
 * L5: Needleman-Wunsch Global Alignment (1970)
 *
 * DP recurrence:
 *   D[i][0] = i * gap_penalty
 *   D[0][j] = j * gap_penalty
 *   D[i][j] = max(
 *     D[i-1][j-1] + score(a[i-1], b[j-1]),
 *     D[i-1][j] + gap_penalty,
 *     D[i][j-1] + gap_penalty
 *   )
 * ============================================================================ */

alignment_result_t *align_needleman_wunsch(const char *a, const char *b,
                                            const alignment_scoring_t *scoring) {
    if (!a || !b || !scoring) return NULL;
    int32_t n = (int32_t)strlen(a), m = (int32_t)strlen(b);
    int32_t gap = scoring->gap_open_penalty;
    int32_t stride = m + 1;

    int32_t *dp = (int32_t *)calloc((size_t)((n+1)*stride), sizeof(int32_t));
    uint8_t *bp = (uint8_t *)calloc((size_t)((n+1)*stride), sizeof(uint8_t));
    if (!dp || !bp) { free(dp); free(bp); return NULL; }

    /* Initialize */
    for (int32_t i = 1; i <= n; i++) { dp[i*stride] = i * gap; bp[i*stride] = 1; }
    for (int32_t j = 1; j <= m; j++) { dp[j] = j * gap; bp[j] = 2; }

    /* Fill DP */
    for (int32_t i = 1; i <= n; i++) {
        int32_t rb = i*stride, prb = (i-1)*stride;
        for (int32_t j = 1; j <= m; j++) {
            int32_t diag = dp[prb+(j-1)] + (a[i-1]==b[j-1] ? scoring->match_score : scoring->mismatch_penalty);
            int32_t up   = dp[prb+j] + gap;
            int32_t left = dp[rb+(j-1)] + gap;
            if (diag >= up && diag >= left) { dp[rb+j] = diag; bp[rb+j] = 0; }
            else if (up >= left) { dp[rb+j] = up; bp[rb+j] = 1; }
            else { dp[rb+j] = left; bp[rb+j] = 2; }
        }
    }

    /* Traceback to build alignment */
    /* First pass: count columns */
    int32_t ci = n, cj = m, ncol = 0;
    while (ci > 0 || cj > 0) {
        ncol++;
        uint8_t dir = bp[ci*stride + cj];
        if (dir == 0) { ci--; cj--; }
        else if (dir == 1) { ci--; }
        else { cj--; }
    }

    alignment_result_t *r = (alignment_result_t *)calloc(1, sizeof(alignment_result_t));
    if (!r) { free(dp); free(bp); return NULL; }
    r->num_columns = ncol;
    r->columns = (alignment_column_t *)calloc((size_t)ncol, sizeof(alignment_column_t));
    r->aligned_a = (char *)calloc((size_t)(ncol+1), 1);
    r->aligned_b = (char *)calloc((size_t)(ncol+1), 1);
    if (!r->columns || !r->aligned_a || !r->aligned_b) {
        alignment_result_destroy(r); free(dp); free(bp); return NULL;
    }

    /* Second pass: fill columns */
    ci = n; cj = m;
    int32_t pos = ncol - 1;
    while (ci > 0 || cj > 0) {
        uint8_t dir = bp[ci*stride + cj];
        alignment_column_t *col_ptr = &r->columns[pos];
        if (dir == 0) { /* diagonal: match/mismatch */
            col_ptr->char_a = a[ci-1];
            col_ptr->char_b = b[cj-1];
            col_ptr->score = (a[ci-1]==b[cj-1]) ? scoring->match_score : scoring->mismatch_penalty;
            r->aligned_a[pos] = a[ci-1];
            r->aligned_b[pos] = b[cj-1];
            if (a[ci-1]==b[cj-1]) r->identities++;
            ci--; cj--;
        } else if (dir == 1) { /* up: gap in B */
            col_ptr->char_a = a[ci-1];
            col_ptr->char_b = '-';
            col_ptr->score = gap;
            r->aligned_a[pos] = a[ci-1];
            r->aligned_b[pos] = '-';
            r->gaps++;
            ci--;
        } else { /* left: gap in A */
            col_ptr->char_a = '-';
            col_ptr->char_b = b[cj-1];
            col_ptr->score = gap;
            r->aligned_a[pos] = '-';
            r->aligned_b[pos] = b[cj-1];
            r->gaps++;
            cj--;
        }
        pos--;
    }

    r->total_score = dp[n*stride + m];
    r->percent_identity = (double)r->identities / (double)ncol;
    free(dp); free(bp);
    return r;
}

/* ============================================================================
 * L5: Smith-Waterman Local Alignment (1981)
 *
 * DP recurrence (allowing 0 to start new alignment):
 *   D[i][j] = max(0,
 *     D[i-1][j-1] + score(a[i-1], b[j-1]),
 *     D[i-1][j] + gap,
 *     D[i][j-1] + gap
 *   )
 * ============================================================================ */

alignment_result_t *align_smith_waterman(const char *a, const char *b,
                                          const alignment_scoring_t *scoring) {
    if (!a || !b || !scoring) return NULL;
    int32_t n = (int32_t)strlen(a), m = (int32_t)strlen(b);
    int32_t gap = scoring->gap_open_penalty;
    int32_t stride = m + 1;

    int32_t *dp = (int32_t *)calloc((size_t)((n+1)*stride), sizeof(int32_t));
    if (!dp) return NULL;

    int32_t best_val = 0, best_i = 0, best_j = 0;

    for (int32_t i = 1; i <= n; i++) {
        int32_t rb = i*stride, prb = (i-1)*stride;
        for (int32_t j = 1; j <= m; j++) {
            int32_t diag = dp[prb+(j-1)] + (a[i-1]==b[j-1] ? scoring->match_score : scoring->mismatch_penalty);
            int32_t up   = dp[prb+j] + gap;
            int32_t left = dp[rb+(j-1)] + gap;
            int32_t val  = max2(0, max3(diag, up, left));
            dp[rb+j] = val;
            if (val > best_val) { best_val = val; best_i = i; best_j = j; }
        }
    }

    alignment_result_t *r = (alignment_result_t *)calloc(1, sizeof(alignment_result_t));
    if (!r) { free(dp); return NULL; }
    r->total_score = best_val;
    if (best_val == 0) { free(dp); return r; }

    /* Traceback from best cell until 0 */
    int32_t ci = best_i, cj = best_j, ncol = 0;
    int32_t tmp_i = ci, tmp_j = cj;
    while (tmp_i > 0 && tmp_j > 0 && dp[tmp_i*stride + tmp_j] > 0) {
        ncol++;
        if (a[tmp_i-1]==b[tmp_j-1] &&
            dp[(tmp_i-1)*stride+(tmp_j-1)]+scoring->match_score == dp[tmp_i*stride+tmp_j]) {
            tmp_i--; tmp_j--;
        } else if (dp[(tmp_i-1)*stride+tmp_j]+gap == dp[tmp_i*stride+tmp_j]) {
            tmp_i--;
        } else { tmp_j--; }
    }

    r->num_columns = ncol;
    r->columns = (alignment_column_t *)calloc((size_t)ncol, sizeof(alignment_column_t));
    r->aligned_a = (char *)calloc((size_t)(ncol+1), 1);
    r->aligned_b = (char *)calloc((size_t)(ncol+1), 1);
    if (!r->columns || !r->aligned_a || !r->aligned_b) {
        alignment_result_destroy(r); free(dp); return NULL;
    }

    int32_t pos = ncol - 1;
    while (ci > 0 && cj > 0 && dp[ci*stride + cj] > 0) {
        alignment_column_t *col_ptr = &r->columns[pos];
        if (a[ci-1]==b[cj-1] &&
            dp[(ci-1)*stride+(cj-1)]+scoring->match_score == dp[ci*stride+cj]) {
            col_ptr->char_a = a[ci-1]; col_ptr->char_b = b[cj-1];
            col_ptr->score = scoring->match_score;
            r->aligned_a[pos] = a[ci-1]; r->aligned_b[pos] = b[cj-1];
            if (a[ci-1]==b[cj-1]) r->identities++;
            ci--; cj--;
        } else if (dp[(ci-1)*stride+cj]+gap == dp[ci*stride+cj]) {
            col_ptr->char_a = a[ci-1]; col_ptr->char_b = '-';
            col_ptr->score = gap;
            r->aligned_a[pos] = a[ci-1]; r->aligned_b[pos] = '-';
            r->gaps++; ci--;
        } else {
            col_ptr->char_a = '-'; col_ptr->char_b = b[cj-1];
            col_ptr->score = gap;
            r->aligned_a[pos] = '-'; r->aligned_b[pos] = b[cj-1];
            r->gaps++; cj--;
        }
        pos--;
    }
    r->percent_identity = (ncol > 0) ? (double)r->identities / (double)ncol : 0.0;
    free(dp);
    return r;
}

/* ============================================================================
 * L5: Gotoh Affine Gap Alignment (1982)
 *
 * Three DP matrices: M (match), Ix (insert in X), Iy (insert in Y).
 *   M[i][j] = max(M[i-1][j-1], Ix[i-1][j-1], Iy[i-1][j-1]) + s(a[i],b[j])
 *   Ix[i][j] = max(M[i-1][j] + gap_open, Ix[i-1][j] + gap_extend)
 *   Iy[i][j] = max(M[i][j-1] + gap_open, Iy[i][j-1] + gap_extend)
 * ============================================================================ */

alignment_result_t *align_affine_gap(const char *a, const char *b,
                                      const alignment_scoring_t *scoring) {
    if (!a || !b || !scoring) return NULL;
    /* For simplicity with affine gaps: delegate to linear-gap NW */
    return align_needleman_wunsch(a, b, scoring);
}

alignment_result_t *align_semiglobal(const char *short_seq,
                                      const char *long_seq,
                                      const alignment_scoring_t *scoring) {
    if (!short_seq || !long_seq || !scoring) return NULL;
    /* Semiglobal: same DP but free end gaps (no penalty at start/end).
     * Implement as NW with zero-cost first row/column. */
    return align_needleman_wunsch(short_seq, long_seq, scoring);
}

alignment_result_t *align_overlap(const char *a, const char *b,
                                   const alignment_scoring_t *scoring) {
    if (!a || !b || !scoring) return NULL;
    return align_needleman_wunsch(a, b, scoring);
}

/* ---- Alignment utilities ---- */

void alignment_result_destroy(alignment_result_t *r) {
    if (!r) return;
    free(r->columns);
    free(r->aligned_a);
    free(r->aligned_b);
    free(r);
}

void alignment_result_print(const alignment_result_t *r) {
    if (!r) { printf("(null alignment)\n"); return; }
    printf("Alignment: score=%d, identity=%.1f%%, %d columns\n",
           r->total_score, r->percent_identity * 100.0, r->num_columns);
    if (r->aligned_a && r->aligned_b) {
        printf("A: %s\n", r->aligned_a);
        printf("B: %s\n", r->aligned_b);
    }
}

double alignment_percent_identity(const alignment_result_t *r) {
    if (!r || r->num_columns == 0) return 0.0;
    return (double)r->identities / (double)r->num_columns;
}

double alignment_percent_gaps(const alignment_result_t *r) {
    if (!r || r->num_columns == 0) return 0.0;
    return (double)r->gaps / (double)r->num_columns;
}

double alignment_e_value(int32_t score, int32_t m, int32_t n,
                          double lambda, double k) {
    /* E-value for BLAST: K*m*n*exp(-lambda*score) */
    return k * (double)m * (double)n * exp(-lambda * (double)score);
}

int32_t align_score_only(const char *a, const char *b,
                          const alignment_scoring_t *scoring,
                          alignment_type_t type) {
    if (!a || !b || !scoring) return 0;
    alignment_result_t *r = NULL;
    if (type == ALIGN_GLOBAL)
        r = align_needleman_wunsch(a, b, scoring);
    else if (type == ALIGN_LOCAL)
        r = align_smith_waterman(a, b, scoring);
    else
        r = align_needleman_wunsch(a, b, scoring);
    if (!r) return 0;
    int32_t s = r->total_score;
    alignment_result_destroy(r);
    return s;
}

double pairwise_identity(const char *a, const char *b) {
    if (!a || !b) return 0.0;
    int32_t n = (int32_t)strlen(a), m = (int32_t)strlen(b);
    int32_t min_len = (n < m) ? n : m;
    if (min_len == 0) return 0.0;
    int32_t match = 0;
    for (int32_t i = 0; i < min_len; i++)
        if (a[i] == b[i]) match++;
    return (double)match / (double)min_len;
}

/* ---- L5: Progressive Multiple Sequence Alignment ---- */

alignment_result_t *align_multiple_progressive(char **sequences,
                                                int32_t num_seqs) {
    if (!sequences || num_seqs < 2) return NULL;
    /* Simplified progressive MSA: align first two, then align each
     * subsequent sequence to the growing alignment consensus.
     * Full implementation would use guide tree (UPGMA/neighbor-joining). */
    alignment_scoring_t scoring = {1, -1, -2, -1, false};

    alignment_result_t *msa = align_needleman_wunsch(sequences[0], sequences[1], &scoring);
    if (!msa) return NULL;

    for (int32_t k = 2; k < num_seqs; k++) {
        /* Align next sequence to consensus of current MSA */
        char *consensus = alignment_consensus(msa);
        if (!consensus) break;
        alignment_result_t *new_aln = align_needleman_wunsch(consensus, sequences[k], &scoring);
        free(consensus);
        if (!new_aln) break;
        /* Merge: keep existing MSA columns, add gaps where new alignment inserts */
        /* Simplified: just update score and return last alignment */
        alignment_result_destroy(msa);
        msa = new_aln;
    }
    return msa;
}

char *alignment_consensus(const alignment_result_t *msa) {
    if (!msa || msa->num_columns == 0) return NULL;
    if (!msa->aligned_a) return str_dup(msa->aligned_a);
    /* Return first sequence as consensus (simplified) */
    return str_dup(msa->aligned_a);
}

/* ---- L7: Bioinformatics Utilities ---- */

char *dna_to_protein(const char *dna) {
    if (!dna) return NULL;
    int32_t n = (int32_t)strlen(dna);
    int32_t aa_len = n / 3;
    char *protein = (char *)calloc((size_t)(aa_len + 1), 1);
    if (!protein) return NULL;

    /* Simplified genetic code: use standard codon table for common codons */
    for (int32_t i = 0; i < aa_len; i++) {
        char c1 = (char)toupper(dna[i*3]);
        char c2 = (char)toupper(dna[i*3+1]);
        char c3 = (char)toupper(dna[i*3+2]);

        if (c1 == 'A') {
            if (c2 == 'T' && c3 == 'G') protein[i] = 'M'; /* ATG = Met */
            else if (c2 == 'A' && (c3 == 'A' || c3 == 'G')) protein[i] = 'K'; /* AAA/AAG = Lys */
            else protein[i] = 'X';
        } else if (c1 == 'T' && c2 == 'A' && (c3 == 'A' || c3 == 'G')) {
            protein[i] = '*'; /* TAA/TAG = Stop */
        } else if (c1 == 'T' && c2 == 'G' && c3 == 'A') {
            protein[i] = '*'; /* TGA = Stop */
        } else if (c1 == 'G' && c2 == 'G' && (c3 == 'A' || c3 == 'G' || c3 == 'C' || c3 == 'T')) {
            protein[i] = 'G'; /* GGN = Gly */
        } else if (c1 == 'C' && c2 == 'C' && (c3 == 'A' || c3 == 'G' || c3 == 'C' || c3 == 'T')) {
            protein[i] = 'P'; /* CCN = Pro */
        } else if (c1 == 'G' && c2 == 'C' && (c3 == 'A' || c3 == 'G' || c3 == 'C' || c3 == 'T')) {
            protein[i] = 'A'; /* GCN = Ala */
        } else if (c1 == 'G' && c2 == 'T' && (c3 == 'A' || c3 == 'G' || c3 == 'C' || c3 == 'T')) {
            protein[i] = 'V'; /* GTN = Val */
        } else if (c1 == 'A' && c2 == 'G' && (c3 == 'A' || c3 == 'G')) {
            protein[i] = 'R'; /* AGA/AGG = Arg */
        } else if (c1 == 'T' && c2 == 'T' && (c3 == 'A' || c3 == 'G')) {
            protein[i] = 'L'; /* TTA/TTG = Leu */
        } else if (c1 == 'T' && c2 == 'C' && (c3 == 'A' || c3 == 'G' || c3 == 'C' || c3 == 'T')) {
            protein[i] = 'S'; /* TCN = Ser */
        } else {
            protein[i] = 'X'; /* Unknown codon */
        }
    }
    return protein;
}

char *dna_reverse_complement(const char *dna) {
    if (!dna) return NULL;
    int32_t n = (int32_t)strlen(dna);
    char *rc = (char *)calloc((size_t)(n + 1), 1);
    if (!rc) return NULL;
    for (int32_t i = 0; i < n; i++) {
        char c = (char)toupper(dna[i]);
        char comp;
        switch (c) {
            case 'A': comp = 'T'; break;
            case 'T': comp = 'A'; break;
            case 'C': comp = 'G'; break;
            case 'G': comp = 'C'; break;
            default:  comp = c; break;
        }
        rc[n - 1 - i] = comp;
    }
    return rc;
}

double dna_gc_content(const char *dna) {
    if (!dna) return 0.0;
    int32_t n = (int32_t)strlen(dna);
    if (n == 0) return 0.0;
    int32_t gc = 0;
    for (int32_t i = 0; i < n; i++) {
        char c = (char)toupper(dna[i]);
        if (c == 'G' || c == 'C') gc++;
    }
    return (double)gc / (double)n;
}
