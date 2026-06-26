#ifndef STRING_UTILS_H
#define STRING_UTILS_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

/* ============================================================================
 * mini-edit-distance-lcs-hardness: String Utilities
 *
 * Helper functions for string manipulation, I/O, and generation
 * used across edit distance and LCS modules.
 *
 * L3: String data structures, suffix/prefix operations
 * L6: Random string generation for benchmarking
 * ============================================================================ */

/* ---- L3: String Utilities ---- */

/* Length of a null-terminated string (like strlen) */
int32_t str_len(const char *s);

/* Copy a string (like strdup) */
char *str_dup(const char *s);

/* Compare two strings for equality */
bool str_eq(const char *a, const char *b);

/* Reverse a string in place */
void str_reverse(char *s);

/* Create reversed copy of a string */
char *str_reversed(const char *s);

/* Check if string is a palindrome */
bool str_is_palindrome(const char *s);

/* Find all occurrences of char in string.
 * Returns array of positions, sets *count. */
int32_t *str_find_all(const char *s, char c, int32_t *count);

/* Compute character frequency histogram.
 * hist[0..255] = count of each ASCII character. */
void str_char_histogram(const char *s, int32_t hist[256]);

/* Compute the alphabet of a string (set of unique chars).
 * Returns null-terminated string of unique chars. */
char *str_alphabet(const char *s);

/* ---- L6: Random String Generation ---- */

/* Generate random string of given length from alphabet.
 * Uses Xorshift PRNG for reproducibility. */
char *str_random(int32_t length, const char *alphabet, uint64_t seed);

/* Generate random DNA sequence (A, C, G, T) */
char *dna_random(int32_t length, uint64_t seed);

/* Generate random protein sequence (20 amino acids) */
char *protein_random(int32_t length, uint64_t seed);

/* Generate random binary string (0, 1) */
char *binary_random(int32_t length, uint64_t seed);

/* Generate a string by mutating an existing string.
 * edit_rate: fraction of characters to edit (0.0 to 1.0).
 * Operations: substitution, insertion, deletion with probabilities. */
char *str_mutate(const char *original, double edit_rate, uint64_t seed);

/* ---- Xorshift PRNG (for reproducibility) ---- */

typedef struct {
    uint64_t state;
} xorshift64_t;

xorshift64_t xorshift64_seed(uint64_t seed);
uint64_t xorshift64_next(xorshift64_t *rng);
int32_t xorshift64_range(xorshift64_t *rng, int32_t min, int32_t max);
double xorshift64_double(xorshift64_t *rng);

#endif /* STRING_UTILS_H */
