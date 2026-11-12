#!/usr/bin/env python3
"""Compare all similarity entries after repeated and changed scoring parameters.

The harness links src/alignment/similarity.c and calls init_sm."""
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

sys.path[:0] = [str(Path(__file__).resolve().parents[1] / name) for name in ('support', 'integration')]
from trf_test_support import SRC, cc

MODULE = 'alignment/similarity.c'
SUPPORT_MODULES = ('core/globals.c', 'core/report.c')


class SimilarityTest(unittest.TestCase):
    def test_cache_and_parameter_changes(self):
        program = r'''
#include <stdlib.h>
#include <limits.h>
#include <assert.h>
#include "alignment/similarity.h"
int main(void) {
    int values[][2]={{0,0},{0,0},{2,-7},{2,-7},{3,-7},{3,-5},
                     {INT_MAX,INT_MIN},{INT_MAX,INT_MIN},{-4,7},{2,-7}};
    int *allocation = NULL;
    for (unsigned n=0;n<sizeof(values)/sizeof(values[0]);++n) {
        int a=values[n][0], b=values[n][1];
        init_sm(a,b);
        if (allocation) assert(allocation==SM);
        allocation=SM;
        for (int i=0;i<256;++i) for (int j=0;j<256;++j) {
            int equal = i==j && (i=='A'||i=='C'||i=='G'||i=='T');
            assert(SM[256*i+j]==(equal?a:b));
            assert(match(i,j)==SM[256*i+j]);
        }
    }
    free(SM);
    return 0;
}
'''
        with tempfile.TemporaryDirectory() as tmp:
            c = Path(tmp) / 'test.c'
            binary = Path(tmp) / 'test'
            c.write_text(program)
            strict = [*cc(), '-O2', '-Wall', '-Wextra', '-Werror', '-I', str(SRC), '-c']
            subprocess.run([*strict, str(c), '-o', str(Path(tmp) / 'test.o')], check=True)
            subprocess.run([*strict, str(SRC / MODULE), '-o', str(Path(tmp) / 'module.o')], check=True)
            subprocess.run([*cc(), '-O2', '-I', str(SRC), str(Path(tmp) / 'test.o'),
                            str(Path(tmp) / 'module.o'),
                            *(str(SRC / name) for name in SUPPORT_MODULES), '-lm',
                            '-o', str(binary)], check=True)
            subprocess.run([str(binary)], check=True)

if __name__ == '__main__':
    unittest.main()
