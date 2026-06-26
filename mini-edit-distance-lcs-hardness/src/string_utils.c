/* ============================================================================
 * string_utils.c -- String Utility Functions
 * L3: String operations  L6: Random string generation
 * ============================================================================ */

#include "string_utils.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

int32_t str_len(const char *s) { return s ? (int32_t)strlen(s) : 0; }

char *str_dup(const char *s) {
    if (!s) return NULL;
    size_t len = strlen(s);
    char *d = (char *)malloc(len + 1);
    if (!d) return NULL;
    memcpy(d, s, len + 1);
    return d;
}

bool str_eq(const char *a, const char *b) {
    if (!a || !b) return (a == b);
    return strcmp(a, b) == 0;
}

void str_reverse(char *s) {
    if (!s) return;
    int32_t n = (int32_t)strlen(s);
    for (int32_t i = 0; i < n / 2; i++) {
        char tmp = s[i]; s[i] = s[n-1-i]; s[n-1-i] = tmp;
    }
}

char *str_reversed(const char *s) {
    if (!s) return NULL;
    char *r = str_dup(s);
    if (r) str_reverse(r);
    return r;
}

bool str_is_palindrome(const char *s) {
    if (!s) return false;
    int32_t n = (int32_t)strlen(s);
    for (int32_t i = 0; i < n/2; i++)
        if (s[i] != s[n-1-i]) return false;
    return true;
}

int32_t *str_find_all(const char *s, char c, int32_t *count) {
    if (!s || !count) return NULL;
    int32_t n = (int32_t)strlen(s), cnt = 0;
    for (int32_t i = 0; i < n; i++) if (s[i] == c) cnt++;
    *count = cnt;
    if (cnt == 0) return NULL;
    int32_t *pos = (int32_t *)malloc((size_t)cnt * sizeof(int32_t));
    if (!pos) return NULL;
    cnt = 0;
    for (int32_t i = 0; i < n; i++) if (s[i] == c) pos[cnt++] = i;
    return pos;
}

void str_char_histogram(const char *s, int32_t hist[256]) {
    if (!s || !hist) return;
    memset(hist, 0, 256 * sizeof(int32_t));
    for (int32_t i = 0; s[i]; i++) hist[(unsigned char)s[i]]++;
}

char *str_alphabet(const char *s) {
    if (!s) return NULL;
    bool seen[256] = {false};
    int32_t count = 0;
    for (int32_t i = 0; s[i]; i++) {
        unsigned char c = (unsigned char)s[i];
        if (!seen[c]) { seen[c] = true; count++; }
    }
    char *alpha = (char *)calloc((size_t)(count + 1), 1);
    if (!alpha) return NULL;
    int32_t idx = 0;
    for (int32_t i = 0; i < 256; i++)
        if (seen[i]) alpha[idx++] = (char)i;
    return alpha;
}

xorshift64_t xorshift64_seed(uint64_t seed) {
    xorshift64_t rng;
    rng.state = (seed == 0) ? 1 : seed;
    return rng;
}

uint64_t xorshift64_next(xorshift64_t *rng) {
    uint64_t x = rng->state;
    x ^= x << 13; x ^= x >> 7; x ^= x << 17;
    rng->state = x;
    return x;
}

int32_t xorshift64_range(xorshift64_t *rng, int32_t min, int32_t max) {
    if (min >= max) return min;
    uint64_t range = (uint64_t)(max - min + 1);
    return min + (int32_t)(xorshift64_next(rng) % range);
}

double xorshift64_double(xorshift64_t *rng) {
    return (double)xorshift64_next(rng) / (double)UINT64_MAX;
}

char *str_random(int32_t length, const char *alphabet, uint64_t seed) {
    if (length < 0 || !alphabet) return NULL;
    int32_t alen = (int32_t)strlen(alphabet);
    if (alen == 0) return NULL;
    char *s = (char *)calloc((size_t)(length + 1), 1);
    if (!s) return NULL;
    xorshift64_t rng = xorshift64_seed(seed);
    for (int32_t i = 0; i < length; i++)
        s[i] = alphabet[(int32_t)(xorshift64_next(&rng) % (uint64_t)alen)];
    return s;
}

char *dna_random(int32_t length, uint64_t seed) {
    return str_random(length, "ACGT", seed);
}

char *protein_random(int32_t length, uint64_t seed) {
    return str_random(length, "ACDEFGHIKLMNPQRSTVWY", seed);
}

char *binary_random(int32_t length, uint64_t seed) {
    return str_random(length, "01", seed);
}

char *str_mutate(const char *original, double edit_rate, uint64_t seed) {
    if (!original || edit_rate < 0.0 || edit_rate > 1.0) return NULL;
    int32_t n = (int32_t)strlen(original);
    xorshift64_t rng = xorshift64_seed(seed);
    int32_t max_len = n + (int32_t)((double)n * edit_rate) + 1;
    char *result = (char *)calloc((size_t)(max_len + 1), 1);
    if (!result) return NULL;
    int32_t ri = 0;
    for (int32_t i = 0; i < n; i++) {
        double r = xorshift64_double(&rng);
        if (r < edit_rate / 3.0) {
            result[ri++] = (char)('a' + (int32_t)(xorshift64_next(&rng) % 26));
        } else if (r < 2.0 * edit_rate / 3.0) {
            result[ri++] = (char)('a' + (int32_t)(xorshift64_next(&rng) % 26));
            result[ri++] = original[i];
        } else if (r < edit_rate) {
            /* deletion */
        } else {
            result[ri++] = original[i];
        }
    }
    result[ri] = 0;
    return result;
}
