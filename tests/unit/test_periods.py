#!/usr/bin/env python3
"""Differential checks against GetTopPeriods in the original TRF source.

Run python3 tests/unit/test_periods.py [--benchmark]. A C compiler is required.
The module src/periods/period_ranking.c is compiled as its own translation
unit and linked with the original function. The trace build adds one
observation call to each, and compares every raw distance count and the
regression slope as well as the returned periods. Allocations of both are
counted through the shim header.
"""

import argparse
import ctypes
from pathlib import Path
import random
import re
import statistics
import subprocess
import sys
import tempfile
import time
import unittest

sys.path[:0] = [str(Path(__file__).resolve().parents[1] / name) for name in ('support', 'integration')]
from trf_test_support import SRC, cc, original_source, require_oracle

MODULE = SRC / "periods/period_ranking.c"
PERIODS = 5
INDEX = [0] * 256
for base, value in zip(b"ACGT", range(4)):
    INDEX[base] = value


SHIM = """/* Counted allocation and histogram observation for both implementations. */
#ifndef PERIODS_SHIM_H
#define PERIODS_SHIM_H
#include <stdlib.h>
#include <string.h>
void *checked_malloc(size_t n);
void *checked_calloc(size_t n, size_t s);
void checked_free(void *p);
void observe(const double *counts, int length, double s);
#define malloc checked_malloc
#define calloc checked_calloc
#define free checked_free
#endif
"""

SUPPORT = """
#include <stdlib.h>
#include <string.h>
int fail_allocation;
static int allocation_count;
static void *live[4];
static size_t live_sizes[4], live_bytes, peak_bytes;
static void *track(void *p, size_t n) {
    int i;
    if (p == NULL) return NULL;
    for (i = 0; i < 4; i++) if (live[i] == NULL) {
        live[i] = p;
        live_sizes[i] = n;
        live_bytes += n;
        if (live_bytes > peak_bytes) peak_bytes = live_bytes;
        return p;
    }
    abort();
}
void *checked_malloc(size_t n) {
    return ++allocation_count == fail_allocation ? NULL : track(malloc(n), n);
}
void *checked_calloc(size_t n, size_t s) {
    return ++allocation_count == fail_allocation ? NULL : track(calloc(n, s), n*s);
}
void checked_free(void *p) {
    int i;
    if (p == NULL) return;
    for (i = 0; i < 4; i++) if (live[i] == p) {
        live[i] = NULL;
        live_bytes -= live_sizes[i];
        free(p);
        return;
    }
    abort();
}
size_t get_live_bytes(void) { return live_bytes; }
size_t get_peak_bytes(void) { return peak_bytes; }
int get_allocation_count(void) { return allocation_count; }
static double *snapshot;
static double slope;
void observe(const double *counts, int length, double s) {
    snapshot = realloc(snapshot, length * sizeof(double));
    memcpy(snapshot, counts, length * sizeof(double));
    slope = s;
}
const double *get_snapshot(void) { return snapshot; }
double get_slope(void) { return slope; }
#include "periods_shim.h"
int *Index;
int MAXDISTANCE;
#define MAXDISTANCECONSTANT 2000
#define NUMBER_OF_PERIODS 5
"""

WRAPPER = """
#include "periods/period_ranking.h"
int reference(const unsigned char *p, int n, int *out, int *index, int max) {
    Index=index; MAXDISTANCE=max; allocation_count=0; peak_bytes=0;
    return GetTopPeriods((unsigned char *)p,n,out);
}
int optimized(const unsigned char *p, int n, int *out, int *index, int max) {
    allocation_count=0; peak_bytes=0;
    return trf_get_top_periods(p,n,out,index,max,5);
}
"""


def build(directory, trace, module_source=None):
    """Shared library of two translation units: the harness with the original
    GetTopPeriods, and the period-ranking module."""
    original = original_source("src/tr30dat.c")
    original = original[original.index("int GetTopPeriods("):]
    original = original[:original.index("\nint multiples_criteria_4(")]
    module = module_source if module_source is not None else MODULE.read_text()
    if trace:
        def instrument(source):
            source, count = re.subn(r"\bs\s*=\s*\(n\s*\*[^;]+;",
                                    lambda m: m[0] + "\nobserve(counts, length, s);",
                                    source)
            assert count == 1
            return source
        original = instrument(original)
        module = instrument(module)
    directory = Path(directory)
    (directory / "periods_shim.h").write_text(SHIM)
    harness = directory / "harness.c"
    harness.write_text(SUPPORT + original + WRAPPER)
    # The module is compiled unchanged apart from the observation call; the
    # shim header, included first, counts its allocations.
    candidate = directory / "period_ranking.c"
    candidate.write_text(module)
    library = directory / "periods.so"
    common = [*cc(), "-O2", "-fPIC", "-I", str(SRC), "-I", str(directory), "-c"]
    subprocess.run([*common, str(harness), "-o", str(directory / "harness.o")], check=True)
    subprocess.run([*common, "-include", "periods_shim.h", str(candidate),
                    "-o", str(directory / "period_ranking.o")], check=True)
    subprocess.run([*cc(), "-shared", str(directory / "harness.o"),
                    str(directory / "period_ranking.o"), "-o", str(library)], check=True)
    lib = ctypes.CDLL(str(library))
    for name in ("reference", "optimized"):
        function = getattr(lib, name)
        function.argtypes = [ctypes.c_char_p, ctypes.c_int,
                             ctypes.POINTER(ctypes.c_int),
                             ctypes.POINTER(ctypes.c_int), ctypes.c_int]
        function.restype = ctypes.c_int
    lib.get_snapshot.restype = ctypes.POINTER(ctypes.c_double)
    lib.get_slope.restype = ctypes.c_double
    lib.get_live_bytes.restype = ctypes.c_size_t
    lib.get_peak_bytes.restype = ctypes.c_size_t
    return lib


def invoke(lib, name, pattern, maxdistance, index=INDEX, trace=True):
    output = (ctypes.c_int * PERIODS)(*([-99] * PERIODS))
    table = (ctypes.c_int * 256)(*index)
    status = getattr(lib, name)(pattern, len(pattern), output, table, maxdistance)
    counts = (ctypes.string_at(lib.get_snapshot(), len(pattern) * 8)
              if trace and status == 0 else None)
    slope = ctypes.c_double(lib.get_slope())
    return status, tuple(output), counts, bytes(slope)


class PeriodsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.directory = tempfile.TemporaryDirectory()
        cls.lib = build(cls.directory.name, True)

    @classmethod
    def tearDownClass(cls):
        cls.directory.cleanup()

    def compare(self, pattern, maxdistance=2000, index=INDEX):
        expected = invoke(self.lib, "reference", pattern, maxdistance, index)
        actual = invoke(self.lib, "optimized", pattern, maxdistance, index)
        self.assertEqual(actual, expected,
                         (len(pattern), maxdistance, pattern[:40]))
        return actual

    def test_short_patterns_and_distance_limits(self):
        for length in range(2, 90):
            for maximum in (0, 1, 2, 7, length - 2, length + 5):
                self.compare((b"ACGTACAA" * 12)[:length], maximum)

    def test_identical_and_equal_scores(self):
        for length in (2, 3, 4, 9, 31, 127, 5999, 6000, 6001, 6002, 6003, 6010):
            result = self.compare(b"A" * length)
            if 9 <= length <= 6002:
                self.assertEqual(result[1], (1, 2, 3, 4, 5))

    def test_cutoff_includes_one_predecessor_at_or_beyond_6000(self):
        for gap in (5999, 6000, 6001, 6007, 12000):
            pattern = b"CC" + b"A" * (gap - 2) + b"CC"
            result = self.compare(pattern, len(pattern))
            counts = (ctypes.c_double * len(pattern)).from_buffer_copy(result[2])
            self.assertEqual(counts[gap], 1.0)
        result = self.compare(b"A" * 6003, 6003)
        counts = (ctypes.c_double * 6003).from_buffer_copy(result[2])
        self.assertEqual(counts[6000], 2.0)
        self.assertEqual(counts[6001], 0.0)
        # Repeated sparse occurrences exercise movement of the cutoff across
        # both exact-boundary and beyond-boundary predecessors.
        for gaps in ((5999, 6000, 6001), (6007, 12000, 5999), (1, 6000, 1)):
            positions = [0]
            for gap in gaps:
                positions.append(positions[-1] + gap)
            pattern = bytearray(b"A" * (positions[-1] + 2))
            for position in positions:
                pattern[position:position + 2] = b"CC"
            self.compare(bytes(pattern), len(pattern))

    def test_exact_integer_histogram_and_peak_storage(self):
        length = 12003
        result = self.compare(b"A" * length, length)
        counts = (ctypes.c_double * length).from_buffer_copy(result[2])
        self.assertEqual(counts[0], 0.0)
        for distance in range(1, length):
            expected = length - 1 - distance if distance <= 6000 else 0
            self.assertEqual(counts[distance], float(expected))
        # The AVX2 histogram replaces the int[length] position buffer with one
        # transient byte buffer (two code copies plus fixed padding/counters);
        # the scalar path keeps the original double[length] + int[length].
        self.assertEqual(self.lib.get_allocation_count(), 2)
        self.assertEqual(self.lib.get_live_bytes(), 0)
        self.assertLessEqual(self.lib.get_peak_bytes(), 12 * length + 32768)

    def test_short_predecessor_ranges_and_cutoff_jump(self):
        # The isolated CC bucket exercises 0 through 7 predecessors before
        # any cutoff applies. A later occurrence advances the lower bound
        # past several predecessors in one step, retaining the nearest one.
        for occurrences in range(1, 9):
            self.compare((b"CCA" * occurrences)[:-1])
        self.compare(b"CCA" * 8 + b"A" * 12000 + b"CC", 20000)

    def test_ambiguous_byte_mapping(self):
        self.compare(b"ANCGTN\x00\xffNACGTN" * 500)
        index = [(i * 7 + 1) % 4 for i in range(256)]
        self.compare(bytes(range(256)) * 30, index=index)

    def test_deterministic_random(self):
        rng = random.Random(355)
        for trial in range(100):
            length = rng.randrange(2, 1600)
            self.compare(bytes(rng.choices(b"ACGTN", k=length)),
                         rng.choice((1, 5, 100, 2000, length + 3)))
        for length in (5998, 6000, 6002, 6200, 12001, 46340, 46343, 50000):
            # Above 46340, deliberately compare the inherited integer i*i
            # expression under the same compiler, without defining new math.
            self.compare(bytes(rng.choices(b"ACGT", k=length)))

    def test_allocation_failure_leaves_results_untouched(self):
        fail = ctypes.c_int.in_dll(self.lib, "fail_allocation")
        try:
            for name, total in (("reference", 2), ("optimized", 2)):
                for allocation in range(1, total + 1):
                    fail.value = allocation
                    result = invoke(self.lib, name, b"ACGT" * 20, 2000)
                    self.assertEqual(result[:2], (1, (-99,) * PERIODS))
                    self.assertEqual(self.lib.get_allocation_count(), allocation)
                    self.assertEqual(self.lib.get_live_bytes(), 0)
                    self.assertLessEqual(self.lib.get_peak_bytes(), 12 * 80)
        finally:
            fail.value = 0


def benchmark():
    rng = random.Random(355)
    patterns = {
        "homopolymer-12000": b"A" * 12000,
        "tandem-12000": (b"ACGTCCATGACTGACATGCA" * 600),
        "random-50000": bytes(rng.choices(b"ACGT", k=50000)),
    }
    with tempfile.TemporaryDirectory() as directory:
        lib = build(directory, False)
        for label, pattern in patterns.items():
            samples = {}
            for name in ("reference", "optimized"):
                times = []
                for _ in range(7):
                    start = time.perf_counter()
                    for repeat in range(3):
                        invoke(lib, name, pattern, 2000, trace=False)
                    times.append((time.perf_counter() - start) / 3)
                samples[name] = statistics.median(times)
            print(f"{label}: original {samples['reference']:.6f}s, "
                  f"NuTRF {samples['optimized']:.6f}s, "
                  f"speedup {samples['reference']/samples['optimized']:.2f}x")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--benchmark", action="store_true")
    options, remaining = parser.parse_known_args()
    require_oracle("source.tar")
    if options.benchmark:
        benchmark()
    else:
        unittest.main(argv=[__file__] + remaining)
