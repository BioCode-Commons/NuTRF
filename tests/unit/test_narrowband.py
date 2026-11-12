#!/usr/bin/env python3
"""Prove the scalar narrowband kernel's score addressing and wraps against the original TRF.

The default unittest suite checks that the kernel source differs from the
original function only by an allowlist of expressions, exhaustively checks
valid columns for periods 1..4000 (including consensus patterns up to twice
MaxPeriod), then compiles the original sources with only narrowbandwrap replaced
by src/alignment/narrowband_scalar.c and compares complete end-to-end output.
A function harness also compares every global state byte on 6000 seeded calls,
including both phases, recentering, and untouched scratch tails.
Use --build-dir PATH to retain the narrowband-only binary for external benchmarks.
Sanitizer findings in newwrap and GetTopPeriods of the original are outside this test.
"""

import argparse
import os
import re
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

sys.path[:0] = [str(Path(__file__).resolve().parents[1] / name) for name in ('support', 'integration')]
from exact_check import compare, execute
from trf_corpus import exact_cases
from trf_test_support import (CFLAGS, SRC, cc, define_guards, extract_original_sources,
                              oracle_dir, original_source, require_oracle, require_sanitizers)
from test_narrowband_state import KERNEL_HEADERS, block_bounds

SCALAR = SRC / 'alignment/narrowband_scalar.c'
SCALAR_HEADER = SRC / 'alignment/narrowband_scalar.h'
DECREMENT = 'c=(c==0)?size-1:c-1;'
INCREMENT = 'c=(c+1==size)?0:c+1;'
REMAINDER_HELPER = """/* Fast C remainder for the common [0, 2*size) row-coordinate range.
 * Negative sentinel expressions and larger values retain C's signed remainder.
 * size is positive; value-size cannot overflow after the nonnegative check. */
static inline int trf_narrowband_remainder(int value, int size)
{
    if(value>=0 && value-size<size)
        return value<size?value:value-size;
    return value%size;
}

"""
SCALAR_INCLUDES = ('#include "core/globals.h"\n#include "core/report.h"\n'
                   '#include "alignment/narrowband_scalar.h"\n')
SCALAR_PROTOTYPE = ('void trf_narrowband_scalar(int start, int size, int bandradius, int bandradiusforward,\n'
                    '                           int option, int tuplesize);\n')


def kernel_text():
    """The remainder helper and row macros of the header followed by the function
    of the source file: the same code, in the same order, as one block."""
    header = SCALAR_HEADER.read_text()
    begin = header.index(REMAINDER_HELPER)
    end = header.index(SCALAR_PROTOTYPE)
    source = SCALAR.read_text()
    assert source.count(SCALAR_INCLUDES) == 1
    function = source[source.index(SCALAR_INCLUDES) + len(SCALAR_INCLUDES):]
    assert function.count('void trf_narrowband_scalar(') == 1
    function = function.replace('void trf_narrowband_scalar(', 'void narrowbandwrap(')
    return header[begin:end].rstrip('\n') + '\n' + function


def build(directory):
    """Compile the original sources with only the narrowband function replaced."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    source_dir = extract_original_sources(directory / 'src')
    source = (source_dir / 'tr30dat.c').read_text()
    start, end = block_bounds(source)
    # The original translation unit already defines the program state, so the
    # module's includes of the shared headers are disabled.
    replacement = (define_guards(KERNEL_HEADERS) +
                   '#define trf_narrowband_scalar narrowbandwrap\n'
                   '#include "alignment/narrowband_scalar.c"\n'
                   '#undef trf_narrowband_scalar\n')
    (source_dir / 'tr30dat.c').write_text(source[:start] + replacement + source[end:])
    binary = directory / 'trf'
    command = [*cc(), *CFLAGS, '-I', str(SRC), str(source_dir / 'trf.c'), '-lm', '-o', str(binary)]
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (directory / 'build.log').write_bytes(result.stdout)
    if result.returncode:
        raise RuntimeError(f'compile failed; see {directory / "build.log"}')
    return binary


class NarrowbandTest(unittest.TestCase):
    def test_kernel_differs_only_by_allowlisted_expressions(self):
        source = original_source('src/tr30dat.c')
        start, end = block_bounds(source)
        header = source[:source.index('*/') + 2]
        self.assertTrue(SCALAR.read_text().startswith(header + '\n\n'))
        candidate = kernel_text().removesuffix('\n')
        self.assertEqual(candidate.count(DECREMENT), 6)
        self.assertEqual(candidate.count(INCREMENT), 6)
        restored = candidate
        self.assertEqual(restored.count(REMAINDER_HELPER), 1)
        restored = restored.replace(REMAINDER_HELPER, '')
        self.assertEqual(restored.count('trf_narrowband_remainder('), 12)
        for expression, count in (
            ('matchatmax_col-1+size', 1), ('Bandcenter[r+1]-1+size', 1),
            ('Bandcenter[r]-Bandcenter[r+1]+size', 1), ('Bandcenter[r]+w', 1),
            ('matchatmax_col-lastmatchatmax_col+size', 2),
            ('mincolbandcenter-1+size', 1), ('matchatmax_col+1', 1),
            ('Bandcenter[r-1]+1', 1), ('Bandcenter[r]-Bandcenter[r-1]+size', 1),
            ('Bandcenter[r]-w+size', 2)):
            call = 'trf_narrowband_remainder(' + expression + ',size)'
            self.assertEqual(restored.count(call), count)
            restored = restored.replace(call, '(' + expression + ')%size')
        comment = ('\t/* Diag/Up originally copy the preceding S row and S+Delta. Each\n'
                   '\t * diagonal index is consumed once, so the read-only row plus local\n'
                   '\t * diag_score preserves the immediate match observer without stores. */\n')
        self.assertEqual(restored.count(comment), 1)
        restored = restored.replace(comment, '')
        declaration = ('\tconst int *previous, *previous_diag, *previous_up;\n'
                       '\tint diag_score;\n')
        self.assertEqual(restored.count(declaration), 1)
        restored = restored.replace(declaration, '')
        self.assertEqual(restored.count('(*pcurr==diag_score)'), 1)
        restored = restored.replace('(*pcurr==diag_score)', '(*pcurr==*pdiag)')
        for indent, offset, comment, following in (
            ('\t\t', 'r+1',
             '\t/* Restore backward scratch once; forward may initialize a narrower band. */\n',
             '\t\tif((matchatmax_col-lastmatchatmax_col+size)%size==size-1)'),
            ('\t\t\t', 'r-1',
             '\t\t/* Restore the caller-visible scratch from the last forward row. */\n',
             '\n\t\t\tif((matchatmax_col-lastmatchatmax_col+size)%size==1)')):
            row_start = restored.index(indent + 'similarity_row=SM+256*currchar;')
            row_end = restored.index(comment, row_start)
            row = restored[row_start:row_end]
            assignment = indent + 'previous=S[' + offset + '];\n'
            self.assertEqual(row.count(assignment), 1)
            row = row.replace(assignment, '')
            self.assertEqual(row.count('diag_score=*previous_diag+(match_yes_no=similarity_row[EC[c]]);'), 3)
            row = row.replace('diag_score=*previous_diag+(match_yes_no=similarity_row[EC[c]]);',
                              '*pdiag+=(match_yes_no=similarity_row[EC[c]]);')
            self.assertEqual(row.count('0,diag_score,'), 3)
            row = row.replace('0,diag_score,', '0,*pdiag,')
            self.assertEqual(row.count('*previous_up+Delta'), 3)
            row = row.replace('*previous_up+Delta', '*pup')
            row = re.sub(r'previous_diag=&previous\[([^\n]+)\];', r'pdiag=&Diag[\1];', row)
            row = re.sub(r'previous_up=&previous\[([^\n]+)\];', r'pup=&Up[\1];', row)
            row = row.replace('previous_diag', 'pdiag').replace('previous_up', 'pup')
            copy = (indent + 'pcurr=&S[r][0];\n' + indent + 'pdiag=&Diag[0];\n' +
                    indent + 'pup=&Up[0];\n' + indent + 'for (i=0;i<=2*w;i++)\n' +
                    indent + '{\n' + indent + '\t*pup=(*pdiag=*pcurr)+Delta;\n' +
                    indent + '\tpcurr++;\n' + indent + '\tpdiag++;\n' +
                    indent + '\tpup++;\n' + indent + '}\n')
            moved_copy = ''.join(line[1:] for line in copy.splitlines(True))
            self.assertEqual(restored.count(comment + moved_copy), 1)
            self.assertEqual(row.count(following), 1)
            row = row.replace(following, copy + following)
            restored = restored[:row_start] + row + restored[row_end:]
            restored = restored.replace(comment + moved_copy, '')
        declaration = '\tconst int *similarity_row;\n'
        self.assertEqual(restored.count(declaration), 1)
        restored = restored.replace(declaration, '')
        for indent in ('\t\t', '\t\t\t'):
            assignment = indent + 'similarity_row=SM+256*currchar;\n'
            self.assertEqual(restored.count('\n' + assignment), 1)
            # Match complete lines so indentation is part of the allowlist.
            restored = restored.replace('\n' + assignment, '\n')
        self.assertEqual(restored.count('similarity_row[EC[c]]'), 6)
        restored = restored.replace('similarity_row[EC[c]]', 'match(currchar,EC[c])')
        restored = restored.replace(DECREMENT, 'c=(c-1+size)%size;')
        restored = restored.replace(INCREMENT, 'c=(c+1)%size;')
        # This byte comparison protects cell order, ties, traceback, score tests,
        # and every recenter/sentinel modulo expression, not just final scores.
        self.assertEqual(restored, source[start:end])

    def test_row_remainder_equivalence(self):
        with tempfile.TemporaryDirectory(prefix='trf-remainder-proof-') as directory:
            source = Path(directory) / 'remainder.c'
            binary = Path(directory) / 'remainder'
            source.write_text('#include <limits.h>\n' + REMAINDER_HELPER + r'''
static int check(int value, int size) {
    return trf_narrowband_remainder(value,size) != value%size;
}
int main(void) {
    for (int size=1;size<=4000;++size) {
        for (int value=-size-3;value<=3*size+3;++value)
            if (check(value,size)) return 1;
        int edges[]={INT_MIN,INT_MIN+1,INT_MAX,INT_MAX-1,-1,0,1};
        for (unsigned int i=0;i<sizeof(edges)/sizeof(edges[0]);++i)
            if (check(edges[i],size)) return 2;
    }
    int divisors[]={INT_MAX,INT_MAX-1,INT_MAX/2,INT_MAX/2+1};
    int values[]={INT_MIN,INT_MIN+1,-1,0,1,INT_MAX/2,INT_MAX-1,INT_MAX};
    for (unsigned int i=0;i<sizeof(divisors)/sizeof(divisors[0]);++i)
        for (unsigned int j=0;j<sizeof(values)/sizeof(values[0]);++j)
            if (check(values[j],divisors[i])) return 3;
    return 0;
}
''')
            subprocess.run([*cc(), '-O1',
                            '-fsanitize=undefined', '-fno-sanitize-recover=all',
                            str(source), '-o', str(binary)], check=True)
            subprocess.run([str(binary)], check=True)

    def test_similarity_row_equivalence(self):
        # A padded matrix also checks every signed-char value without forming
        # an out-of-bounds pointer. TRF's normalized input uses positive ASCII.
        with tempfile.TemporaryDirectory(prefix='trf-score-row-proof-') as directory:
            source = Path(directory) / 'row.c'
            binary = Path(directory) / 'row'
            source.write_text(r'''#include <limits.h>
int main(void) {
    int storage[131072];
    int *SM = storage + 65536;
    for (int value=CHAR_MIN; value<=CHAR_MAX; ++value) {
        char currchar=value;
        const int *similarity_row=SM+256*currchar;
        for (int ec=0; ec<=UCHAR_MAX; ++ec) {
            if (&similarity_row[ec] != &SM[256*currchar+ec]) return 1;
        }
    }
    return 0;
}
''')
            subprocess.run([*cc(), '-O0',
                            '-fsanitize=undefined', '-fno-sanitize-recover=all',
                            str(source), '-o', str(binary)], check=True)
            subprocess.run([str(binary)], check=True)

    def test_exhaustive_wrap_equivalence(self):
        # Execute C expressions so Python's different negative-modulo convention
        # cannot hide an issue. Only the canonical 0 <= c < size domain is valid.
        with tempfile.TemporaryDirectory(prefix='trf-wrap-proof-') as directory:
            source = Path(directory) / 'wrap.c'
            binary = Path(directory) / 'wrap'
            source.write_text('''
int main(void) {
    for (int size=1; size<=4000; ++size) {
        for (int initial=0; initial<size; ++initial) {
            int c=initial;
            int old_down=(c-1+size)%size;
            ''' + DECREMENT + '''
            if(c!=old_down || c<0 || c>=size) return 1;
            c=initial;
            int old_up=(c+1)%size;
            ''' + INCREMENT + '''
            if(c!=old_up || c<0 || c>=size) return 2;
        }
    }
    return 0;
}
''')
            subprocess.run([*cc(), '-O0',
                            '-fsanitize=undefined', '-fno-sanitize-recover=all',
                            str(source), '-o', str(binary)], check=True)
            subprocess.run([str(binary)], check=True)

    def test_complete_global_state(self):
        from test_narrowband_state import run
        with tempfile.TemporaryDirectory(prefix='trf-narrowband-state-') as directory:
            run(directory)

    def test_narrowband_only_end_to_end(self):
        names = {'small_test-all', 'small_test-ngs', 'diff_cn-all',
                 'overlap', 'pm75-pi20', 'pm75-pi20-mixed-baseline-failure',
                 'maxperiod-2000', 'alternate-scoring', 'lowercase-ambiguous',
                 'multiple-records', 'pagination-html', 'crlf', 'no-final-newline',
                 'substitution-indel-2', 'substitution-indel-3',
                 'substitution-indel-5'}
        names.update('period-' + str(p) for p in
                     (1, 2, 20, 21, 50, 51, 100, 101, 499, 500, 501, 1999, 2000))
        cases = [case for case in exact_cases() if case.name in names]
        self.assertEqual(len(cases), len(names))
        with tempfile.TemporaryDirectory(prefix='trf-narrowband-') as directory:
            directory = Path(directory)
            binary = build(directory / 'build')
            for case in cases:
                with self.subTest(case=case.name):
                    reference_work = directory / case.name / 'reference'
                    candidate_work = directory / case.name / 'candidate'
                    reference = execute(oracle_dir() / 'trf', case, reference_work, 120)
                    candidate = execute(binary, case, candidate_work, 120)
                    expected_exit = case.expected_exit
                    self.assertEqual(reference['exit'], expected_exit)
                    self.assertEqual(candidate['exit'], expected_exit)
                    self.assertEqual(compare(reference, candidate,
                                             reference_work, candidate_work), [])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build-dir', type=Path,
                        help='retain the narrowband-only binary built from the original sources')
    parser.add_argument('--oracle', type=Path,
                        help='directory containing source.tar and trf of the original; overrides TRF_ORACLE_DIR')
    parser.add_argument('--wrap-only', action='store_true',
                        help='C proofs of the wrap and remainder expressions; needs no original TRF')
    options, remaining = parser.parse_known_args()
    if options.oracle:
        os.environ['TRF_ORACLE_DIR'] = str(options.oracle.resolve())
    require_sanitizers('-fsanitize=undefined')
    if options.wrap_only:
        selected = ['NarrowbandTest.' + name for name in (
            'test_exhaustive_wrap_equivalence', 'test_row_remainder_equivalence',
            'test_similarity_row_equivalence')]
        unittest.main(argv=[__file__] + selected + remaining)
    require_oracle()
    if options.build_dir:
        print(build(options.build_dir.resolve()))
    else:
        unittest.main(argv=[__file__] + remaining)
