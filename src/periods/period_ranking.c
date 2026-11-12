/*
Tandem Repeats Finder
Copyright (C) 1999-2020 Gary Benson

This file is part of the Tandem Repeats Finder (TRF) program.

TRF is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of
the License, or (at your option) any later version.

TRF is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See
the GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public
License along with TRF. If not, see <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>

#include "periods/period_ranking.h"

/* Distances at or beyond this cutoff are counted only for the nearest such
 * predecessor (MAXDISTANCECONSTANT*3 in the original source). */
#define TRF_PERIODS_CUTOFF 6000

#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__)) && \
    !defined(TRF_PERIODS_DISABLE_SIMD)
#define TRF_PERIODS_HAVE_AVX2 1
#include <immintrin.h>
#include <stdint.h>
#include <string.h>

/* Patterns shorter than this keep the scalar pair loop. Both paths produce
 * the same integer histogram, so the threshold only affects speed. */
#ifndef TRF_PERIODS_SIMD_MIN_LENGTH
#define TRF_PERIODS_SIMD_MIN_LENGTH 128
#endif
#define TRF_PERIODS_BLOCK 4096          /* multiple of 32, at most 255*32 */
#define TRF_PERIODS_FRONT (TRF_PERIODS_CUTOFF + 32)
#define TRF_PERIODS_TAIL 64

/* near[d] += #{ i in [d,n) : a[i]==b[i-d] } for 1 <= d < min(n,CUTOFF).
 * a and b hold the same n tuple codes (0..15). a is followed by 0xFF bytes,
 * b is surrounded by 0xFE bytes, so whole 32-byte compares may run past
 * either end of the valid range without ever producing a spurious match.
 * Positions are blocked so that the block and its 6000-byte look-back stay in
 * the first-level cache; four distances share each load of a. A byte counter
 * sees at most BLOCK/32 <= 255 increments before it is widened. */
__attribute__((target("avx2")))
static void trf_periods_near_avx2(const unsigned char *a_codes,
                                  const unsigned char *b_codes, int n,
                                  uint32_t *near)
{
    const __m256i zero = _mm256_setzero_si256();
    int s, e, len, dmax, d, k;

    for (s = 0; s < n; s = e) {
        e = (n - s > TRF_PERIODS_BLOCK) ? s + TRF_PERIODS_BLOCK : n;
        len = (e - s + 31) & ~31;
        dmax = (e - 1 < TRF_PERIODS_CUTOFF - 1) ? e - 1 : TRF_PERIODS_CUTOFF - 1;
        const unsigned char *a = a_codes + s;
        for (d = 1; d + 3 <= dmax; d += 4) {
            const unsigned char *b = b_codes + s - d;
            __m256i c0 = zero, c1 = zero, c2 = zero, c3 = zero, sums;
            __m128i fold;
            uint64_t packed;
            for (k = 0; k < len; k += 32) {
                __m256i va = _mm256_loadu_si256((const __m256i *)(a + k));
                c0 = _mm256_sub_epi8(c0, _mm256_cmpeq_epi8(va,
                        _mm256_loadu_si256((const __m256i *)(b + k))));
                c1 = _mm256_sub_epi8(c1, _mm256_cmpeq_epi8(va,
                        _mm256_loadu_si256((const __m256i *)(b + k - 1))));
                c2 = _mm256_sub_epi8(c2, _mm256_cmpeq_epi8(va,
                        _mm256_loadu_si256((const __m256i *)(b + k - 2))));
                c3 = _mm256_sub_epi8(c3, _mm256_cmpeq_epi8(va,
                        _mm256_loadu_si256((const __m256i *)(b + k - 3))));
            }
            /* Each total is at most BLOCK+31 < 65536: four 16-bit fields. */
            sums = _mm256_sad_epu8(c0, zero);
            sums = _mm256_or_si256(sums, _mm256_slli_epi64(_mm256_sad_epu8(c1, zero), 16));
            sums = _mm256_or_si256(sums, _mm256_slli_epi64(_mm256_sad_epu8(c2, zero), 32));
            sums = _mm256_or_si256(sums, _mm256_slli_epi64(_mm256_sad_epu8(c3, zero), 48));
            fold = _mm_add_epi64(_mm256_castsi256_si128(sums),
                                 _mm256_extracti128_si256(sums, 1));
            fold = _mm_add_epi64(fold, _mm_unpackhi_epi64(fold, fold));
            packed = (uint64_t)_mm_cvtsi128_si64(fold);
            near[d] += (uint32_t)(packed & 0xFFFF);
            near[d + 1] += (uint32_t)((packed >> 16) & 0xFFFF);
            near[d + 2] += (uint32_t)((packed >> 32) & 0xFFFF);
            near[d + 3] += (uint32_t)(packed >> 48);
        }
        for (; d <= dmax; d++) {
            const unsigned char *b = b_codes + s - d;
            __m256i c0 = zero, sums;
            __m128i fold;
            for (k = 0; k < len; k += 32) {
                __m256i va = _mm256_loadu_si256((const __m256i *)(a + k));
                c0 = _mm256_sub_epi8(c0, _mm256_cmpeq_epi8(va,
                        _mm256_loadu_si256((const __m256i *)(b + k))));
            }
            sums = _mm256_sad_epu8(c0, zero);
            fold = _mm_add_epi64(_mm256_castsi256_si128(sums),
                                 _mm256_extracti128_si256(sums, 1));
            fold = _mm_add_epi64(fold, _mm_unpackhi_epi64(fold, fold));
            near[d] += (uint32_t)_mm_cvtsi128_si64(fold);
        }
    }
}

/* Fill counts[] with exactly the histogram of the scalar pair loop below.
 * Returns 0 when unavailable (old CPU, short pattern, allocation failure);
 * counts[] is then untouched. */
static int trf_periods_histogram_simd(const unsigned char *pattern, int length,
                                      const int *index, double *counts)
{
    unsigned char *buffer, *a_codes, *b_codes;
    uint32_t *near;
    int last[16];
    int n = length - 1;         /* tuple positions 0..length-2 */
    int i, d, tupid, dmax;
    size_t a_size, b_size;

    if (length < TRF_PERIODS_SIMD_MIN_LENGTH || length > INT32_MAX - 4096) return 0;
#ifdef TRF_PERIODS_FORCE_SCALAR
    return 0;
#endif
    if (!__builtin_cpu_supports("avx2")) return 0;

    a_size = (size_t)n + TRF_PERIODS_TAIL;
    b_size = (size_t)TRF_PERIODS_FRONT + (size_t)n + TRF_PERIODS_TAIL;
    buffer = (unsigned char *)malloc(a_size + b_size +
                                     TRF_PERIODS_CUTOFF * sizeof(uint32_t));
    if (buffer == NULL) return 0;
    near = (uint32_t *)buffer;      /* malloc alignment suffices */
    a_codes = buffer + TRF_PERIODS_CUTOFF * sizeof(uint32_t);
    b_codes = a_codes + a_size + TRF_PERIODS_FRONT;

    for (i = 0; i < n; i++) {
        tupid = index[pattern[i]] * 4 + index[pattern[i + 1]];
        if (tupid < 0 || tupid > 15) {  /* not the production index table */
            free(buffer);
            return 0;
        }
        a_codes[i] = (unsigned char)tupid;
    }
    memset(a_codes + n, 0xFF, TRF_PERIODS_TAIL);
    memset(b_codes - TRF_PERIODS_FRONT, 0xFE, TRF_PERIODS_FRONT);
    memcpy(b_codes, a_codes, (size_t)n);
    memset(b_codes + n, 0xFE, TRF_PERIODS_TAIL);
    memset(near, 0, TRF_PERIODS_CUTOFF * sizeof(uint32_t));

    /* every predecessor nearer than the cutoff */
    trf_periods_near_avx2(a_codes, b_codes, n, near);
    dmax = (n - 1 < TRF_PERIODS_CUTOFF - 1) ? n - 1 : TRF_PERIODS_CUTOFF - 1;
    for (d = 1; d <= dmax; d++) counts[d] += (double)near[d];

    /* plus the nearest predecessor at or beyond the cutoff, if any */
    for (tupid = 0; tupid < 16; tupid++) last[tupid] = -1;
    for (i = TRF_PERIODS_CUTOFF; i < n; i++) {
        last[a_codes[i - TRF_PERIODS_CUTOFF]] = i - TRF_PERIODS_CUTOFF;
        d = last[a_codes[i]];
        if (d >= 0) counts[i - d] += 1.0;
    }

    free(buffer);
    return 1;
}
#else
#define TRF_PERIODS_HAVE_AVX2 0
#endif

int trf_get_top_periods(const unsigned char *pattern, int length, int *toparray,
                        const int *index, int maxdistance,
                        int number_of_periods)
{
    int sizes[16] = {0};
    int starts[17];
    int next[16];
    int *positions;
    double *counts;
    int i, t, end, tupid, current, previous, lower, dist;
    int topind;
    double topval;
    double n, xysum, xsum, ysum, x2sum, s;

    counts = (double *)calloc(length, sizeof(double));
    if (counts == NULL) return 1;
    positions = NULL;

#if TRF_PERIODS_HAVE_AVX2
    if (trf_periods_histogram_simd(pattern, length, index, counts))
        goto histogram_done;
#endif

    positions = (int *)malloc(length * sizeof(int));
    if (positions == NULL) {
        free(counts);
        return 1;
    }

    /* Each tuple's occurrences occupy one contiguous, ascending bucket. */
    end = length - 2;
    for (i = 0; i <= end; i++) {
        tupid = index[pattern[i]] * 4 + index[pattern[i + 1]];
        sizes[tupid]++;
    }
    starts[0] = 0;
    for (tupid = 0; tupid < 16; tupid++) {
        starts[tupid + 1] = starts[tupid] + sizes[tupid];
        next[tupid] = starts[tupid];
    }
    for (i = 0; i <= end; i++) {
        tupid = index[pattern[i]] * 4 + index[pattern[i + 1]];
        positions[next[tupid]++] = i;
    }

    for (tupid = 0; tupid < 16; tupid++) {
        lower = starts[tupid];
        for (current = starts[tupid] + 1;
             current < starts[tupid + 1]; current++) {
            i = positions[current];
            /* Keep the nearest predecessor at or beyond the cutoff, plus
             * every nearer predecessor. The original includes that first
             * distance >= 6000 even when a sparse gap makes it much larger.
             * Ascending positions make this lower bound move only forward.
             */
            while (lower + 1 < current &&
                   i - positions[lower + 1] >= TRF_PERIODS_CUTOFF)
                lower++;
            for (previous = lower; previous < current; previous++) {
                dist = i - positions[previous];
                counts[dist] += 1.0;
            }
        }
    }

#if TRF_PERIODS_HAVE_AVX2
histogram_done:
#endif
    /* Keep the original order and expressions: in particular, i*i is integer
     * arithmetic. Signed overflow for long patterns is a known upstream defect.
     * Reordering +1 updates above is exact: each distance receives fewer than
     * length increments, well below the exact integer range of double.
     */
    xysum = xsum = ysum = x2sum = 0.0;
    end = length - 2;
    for (i = 1; i <= end; i++) {
        xysum += (i * counts[i]);
        xsum += (i);
        ysum += (counts[i]);
        x2sum += (i * i);
    }
    n = end;
    s = (n * xysum - xsum * ysum) / (n * x2sum - xsum * xsum);

    end = length - 2;
    for (i = 1; i <= end; i++) {
        counts[i] = counts[i] - i * s;
    }

    end = length - 2;
    if (end > maxdistance) end = maxdistance;
    for (t = 0; t < number_of_periods; t++) {
        topind = 0;
        topval = 0.0;
        for (i = 1; i <= end; i++) {
            if (counts[i] > topval) {
                topind = i;
                topval = counts[i];
            }
        }
        toparray[t] = topind;
        counts[topind] = 0.0;
    }

    free(counts);
    free(positions);
    return 0;
}
