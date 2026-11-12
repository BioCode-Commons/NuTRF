#!/usr/bin/env python3
"""Compare repeat-list sort node identity and payloads against the original TRF source.

The harness links src/redundancy/list_sort.c and compiles the original
SortByIndex and SortByCount next to it under other names.
Run python3 tests/unit/test_list_sort.py [--benchmark]. The bounded benchmark
measures list-sort scalability only; it makes no end-to-end TRF speed claim.
"""

import argparse
import ctypes
import os
from pathlib import Path
import random
import re
import shlex
import statistics
import subprocess
import sys
import tempfile
import unittest

sys.path[:0] = [str(Path(__file__).resolve().parents[1] / name) for name in ('support', 'integration')]
from trf_test_support import SRC, cc, original_source, require_oracle

MODULE = SRC / 'redundancy/list_sort.c'
INT_MIN = -(2**31)
INT_MAX = 2**31 - 1


def function(source, name):
    match = re.search(r"(?:static\s+)?IL\*\s+" + name + r"\([^;]*?\)\s*\{", source)
    assert match is not None, name
    position = source.index("{", match.start())
    depth = 1
    end = position + 1
    while depth:
        depth += (source[end] == "{") - (source[end] == "}")
        end += 1
    return source[match.start():end]


def build(directory):
    original = original_source('src/trfclean.h')
    code = "#include <stdlib.h>\n#include <stddef.h>\n#include <stdint.h>\n"
    code += "#include <string.h>\n#include <time.h>\n"
    code += '#include "redundancy/list_sort.h"\n'
    for name in ("SortByIndex", "SortByCount"):
        code += function(original, name).replace(name, "original_" + name)
    code += r"""
static IL *sort(int candidate, int by_count, IL *head) {
    if (candidate) return by_count ? SortByCount(head) : SortByIndex(head);
    return by_count ? original_SortByCount(head) : original_SortByIndex(head);
}

/* Check addresses before dereferencing a returned link, then bound traversal.
 * Compare every payload byte, including padding, excluding only next.
 */
int run_sort(int candidate, int by_count, int n,
             const int *first, const int *count, int *order) {
    IL *nodes, *before, *node;
    unsigned char *seen;
    int i, status = 0;
    uintptr_t address, base, offset;
    size_t id;
    if (n < 0) return 7;
    if (n == 0) return sort(candidate, by_count, NULL) != NULL;
    nodes = malloc(n * sizeof(*nodes));
    before = malloc(n * sizeof(*before));
    seen = calloc(n, 1);
    if (!nodes || !before || !seen) abort();
    for (i = 0; i < n; i++) {
        memset(nodes + i, (i % 251) + 1, sizeof(*nodes));
        nodes[i].first = first[i];
        nodes[i].count = count[i];
        nodes[i].pattern = (char *)(nodes + i);
        nodes[i].next = i + 1 < n ? nodes + i + 1 : NULL;
    }
    memcpy(before, nodes, n * sizeof(*nodes));
    node = sort(candidate, by_count, nodes);
    base = (uintptr_t)nodes;
    for (i = 0; i < n; i++) {
        address = (uintptr_t)node;
        if (address < base || address - base >= n * sizeof(*nodes)) {
            status = 2; break;
        }
        offset = address - base;
        if (offset % sizeof(*nodes) != 0) { status = 3; break; }
        id = offset / sizeof(*nodes);
        if (seen[id]) { status = 4; break; }
        seen[id] = 1;
        if (memcmp(nodes + id, before + id, offsetof(IL, next))) {
            status = 5; break;
        }
        order[i] = id;
        node = node->next;
    }
    if (!status && node != NULL) status = 6;
    free(nodes); free(before); free(seen);
    return status;
}

/* Rebuild links outside the measured section; time only C sorting CPU time. */
double benchmark(int candidate, int by_count, int n, const int *keys,
                 int repeats) {
    IL *nodes = calloc(n, sizeof(*nodes));
    IL *head;
    int i, r;
    clock_t start;
    double elapsed = 0.0;
    volatile int sink = 0;
    if (!nodes) abort();
    for (r = 0; r < repeats; r++) {
        for (i = 0; i < n; i++) {
            nodes[i].first = nodes[i].count = keys[i];
            nodes[i].next = i + 1 < n ? nodes + i + 1 : NULL;
        }
        start = clock();
        head = sort(candidate, by_count, nodes);
        elapsed += (double)(clock() - start) / CLOCKS_PER_SEC;
        sink ^= head->first;
    }
    free(nodes);
    return elapsed / repeats + 0.0 * sink;
}
"""
    source = Path(directory) / "list_sort.c"
    library = Path(directory) / "list_sort.so"
    source.write_text(code)
    subprocess.run(cc() +
                   shlex.split(os.environ.get("LIST_SORT_CFLAGS", "-O2")) +
                   ["-Wall", "-Wextra", "-shared", "-fPIC", "-I", str(SRC),
                    str(source), str(MODULE), "-o", str(library)], check=True)
    lib = ctypes.CDLL(str(library))
    pointer = ctypes.POINTER(ctypes.c_int)
    lib.run_sort.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int,
                            pointer, pointer, pointer]
    lib.run_sort.restype = ctypes.c_int
    lib.benchmark.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int,
                             pointer, ctypes.c_int]
    lib.benchmark.restype = ctypes.c_double
    return lib


class ListSortTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.directory = tempfile.TemporaryDirectory()
        cls.lib = build(cls.directory.name)

    @classmethod
    def tearDownClass(cls):
        cls.directory.cleanup()

    def compare(self, first, count=None):
        if count is None:
            # Independent values prove the comparator ignores the other key.
            count = [(i * 73 + 11) % 31 for i in range(len(first))]
        n = len(first)
        first_array = (ctypes.c_int * n)(*first)
        count_array = (ctypes.c_int * n)(*count)
        for by_count, keys in ((0, first), (1, count)):
            expected = sorted(range(n), key=keys.__getitem__)
            for candidate in (0, 1):
                order = (ctypes.c_int * n)()
                status = self.lib.run_sort(candidate, by_count, n,
                                           first_array, count_array, order)
                self.assertEqual(status, 0, (candidate, by_count, n))
                self.assertEqual(list(order), expected,
                                 (candidate, by_count, n))

    def test_empty_and_singleton(self):
        self.compare([])
        self.compare([INT_MIN], [INT_MAX])

    def test_sorted_reversed_equal_and_alternating(self):
        for n in (2, 3, 4, 7, 8, 9, 31, 32, 33, 127, 128, 129, 2000):
            for keys in (list(range(n)), list(range(n - 1, -1, -1)),
                         [7] * n, [i % 2 for i in range(n)]):
                self.compare(keys)
                self.compare(list(reversed(keys)), keys)

    def test_extreme_keys(self):
        self.compare([INT_MAX, INT_MIN, 0, INT_MIN, INT_MAX, -1, 1] * 21)

    def test_seeded_random_and_ties(self):
        rng = random.Random(355)
        for _ in range(150):
            n = rng.randrange(0, 600)
            self.compare([rng.randint(-15, 15) for _ in range(n)],
                         [rng.randint(INT_MIN, INT_MAX) for _ in range(n)])


def benchmark():
    rng = random.Random(355)
    with tempfile.TemporaryDirectory() as directory:
        lib = build(directory)
        for n in (2000, 10000):
            patterns = {"sorted": list(range(n)),
                        "reversed": list(range(n - 1, -1, -1)),
                        "equal": [7] * n,
                        "random": [rng.randrange(n // 4) for _ in range(n)]}
            for label, keys in patterns.items():
                array = (ctypes.c_int * n)(*keys)
                times = {0: [], 1: []}
                for pair in range(5):
                    for candidate in ((0, 1) if pair % 2 == 0 else (1, 0)):
                        repeats = 500 if label in ("sorted", "equal") else 1
                        times[candidate].append(lib.benchmark(
                            candidate, 0, n, array, repeats))
                old, new = (statistics.median(times[key]) for key in (0, 1))
                print(f"n={n} {label}: original={old:.8f}s "
                      f"merge={new:.8f}s speedup={old/new:.2f}x")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--benchmark", action="store_true")
    options, remaining = parser.parse_known_args()
    require_oracle('source.tar')
    if options.benchmark:
        benchmark()
    else:
        unittest.main(argv=[__file__] + remaining)
