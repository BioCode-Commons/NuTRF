#!/usr/bin/env python3
"""Verify exact newwrap cell reuse and the traceback projection.

The harness holds the program state in small arrays and compiles against it the
original newwrap, the module src/alignment/newwrap.c, and the traceback module
src/alignment/pair_alignment.c, with their includes of the shared headers
disabled. The test-only full-state implementation compares every row of all
three scans, including overwritten rows, sentinels, Up/Diag, maxima and ties.
The skip-middle implementation used by the program compares endpoints, live
traceback matrix, Up/Diag, EC and actual alignment strings/indices/copy counts.
Degenerate cases initialize the original's otherwise indeterminate endpoint
locals identically in both test functions. Periods 1..25 exercise the <=20
specialization and its fallback.

--build-dir retains a binary built from the original sources with only newwrap replaced.
--kernel-only uses the copy of the original function kept in newwrap.c in place
of the original sources, so it needs no original TRF.
TRF_KERNEL_CFLAGS may enable sanitizers; TRF_SANITIZER_PAD_ORACLE adds a leading
buffer element solely to accommodate the original's one-before pointer formation.
"""
import argparse
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys
import tempfile
import unittest

sys.path[:0] = [str(Path(__file__).resolve().parents[1] / name) for name in ('support', 'integration')]
from exact_check import compare, execute
from trf_corpus import exact_cases
from trf_test_support import (CFLAGS, SRC, cc, define_guards, extract_original_sources,
                              oracle_dir, original_source, require_oracle)

MODULE = SRC / 'alignment/newwrap.c'
TRACEBACK = SRC / 'alignment/pair_alignment.c'


def block_bounds(source):
    start = source.index('/***************************** newwrap()')
    start = source.rfind('/*******************************************************************/', 0, start)
    end = source.index('\n}\n', start) + 2
    return start, end


def original_block():
    source = original_source('src/tr30dat.c')
    start, end = block_bounds(source)
    return source[start:end]


def fallback_block():
    source = MODULE.read_text()
    start = source.index('/*******************************************************************/')
    end = source.index('\n}\n', start) + 2
    return source[start:end].replace('static void newwrap_original(', 'void newwrap(', 1)


def build(directory):
    """Compile the original sources with only newwrap replaced by the module."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    source_dir = extract_original_sources(directory / 'src')
    path = source_dir / 'tr30dat.c'
    source = path.read_text()
    start, end = block_bounds(source)
    # The original translation unit already defines the program state, so the
    # module's includes of the shared headers are disabled.
    path.write_text(source[:start] + define_guards() + '#include "alignment/newwrap.c"\n' + source[end:])
    binary = directory / 'trf'
    result = subprocess.run([*cc(), *CFLAGS, '-I', str(SRC),
                             str(source_dir / 'trf.c'), '-lm', '-o', str(binary)],
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (directory / 'build.log').write_bytes(result.stdout)
    if result.returncode:
        raise RuntimeError(f'compile failed; see {directory / "build.log"}')
    return binary


def traced(source):
    source += "\n"
    # These variables are uninitialized in the original's no-positive-match cases.
    # Set only in this harness to make the otherwise undefined cases comparable.
    source = source.replace('maxscore,minrow,maxrow,maxcol,modstart,maxrealrow;',
                            'maxscore,minrow=start,maxrow=0,maxcol=0,modstart,maxrealrow=0;')
    # Instrument each function independently; the projection variant has two scans.
    functions = list(re.finditer(r'(?:static )?void (?:original_newwrap|newwrap\w*)\(', source))
    for function in reversed(functions):
        end = source.index('\n}\n', function.start()) + 2
        body = source[function.start():end]
        matches = list(re.finditer(r'Diag\[(?:adjlength|0)\]=(\*pcurr|S\[r\]\[0\]|backward_row\[0\]|Diag\[size\]);', body))
        assert len(matches) in (2, 3), len(matches)
        phases = (0, 1, 2) if len(matches) == 3 else (0, 2)
        for phase, match in reversed(list(zip(phases, matches))):
            statement = ('\nrecord_row(%d,r,realr,end_of_trace,maxscore,minrow,maxrow,'
                         'maxcol,maxrealrow,pleft);' % phase)
            body = body[:match.end()] + statement + body[match.end():]
        source = source[:function.start()] + body + source[end:]

    return source


SUPPORT = r'''
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#define COLS 26
#define ROWS 99
#define FALSE 0
#define TRUE 1
#define max(a,b) ((a)>=(b)?(a):(b))
#define max2(a,b) max(a,b)
#define max4(a,b,c,d) max(max(a,b),max(c,d))
#define match(a,b) SM[256*(a)+(b)]
static unsigned char sequence_storage[270], *Sequence=sequence_storage, EC[COLS];
static struct {char pattern[COLS];} Consensus;
static int matrix[ROWS][COLS], *rowptrs[ROWS], **S=rowptrs;
#ifdef TRF_SANITIZER_PAD_ORACLE
static int up_storage[COLS+2],diag_storage[COLS+2];
#define Up (up_storage+1)
#define Diag (diag_storage+1)
#else
static int Up[COLS], Diag[COLS];
#endif
static int sm_storage[65536], *SM=sm_storage;
static int Length,Alpha,Beta,Delta,Min_Distance_Window,Criteria_print=0;
static unsigned long avx2_calls;
static unsigned int maxwraplength;
static double Rows;
static int Wrapend,Reportmin,Maxrealrow,Maxrow,Maxcol,Maxscore;
static void trf_message(const char *fmt,...) {(void)fmt;abort();}
struct row_state {int local[10], row[COLS], up[COLS], diag[COLS];};
static struct row_state trace[3*ROWS];
static unsigned int ntrace, trace_index, candidate, testcase, phase_seen[3], traceback_cases;
static int projection;
static void record_row(int phase,unsigned int r,int realr,int end,int score,
                       int minrow,int maxrow,int maxcol,int maxreal,int left) {
    if(projection) return;
    struct row_state now={{phase,(int)r,realr,end,score,minrow,maxrow,maxcol,maxreal,left},{0},{0},{0}};
    memcpy(now.row,S[r],sizeof now.row);memcpy(now.up,Up,sizeof now.up);memcpy(now.diag,Diag,sizeof now.diag);
    if(!candidate) {if(ntrace>=3*ROWS) abort();trace[ntrace++]=now;phase_seen[phase]++;}
    else if(trace_index>=ntrace || memcmp(&now,&trace[trace_index],sizeof now)) {
        fprintf(stderr,"case %u phase %d row %u trace %u differs (size parameter in case log)\n",testcase,phase,r,trace_index);exit(1);
    }
    trace_index++;
}
static unsigned int rngstate;
static unsigned int rnd(void) {rngstate^=rngstate<<13;rngstate^=rngstate>>17;return rngstate^=rngstate<<5;}
static void setup(unsigned int seed,int size,int style,int start,int length,int wrap,int matchscore,int mismatch,int gap,int window) {
    rngstate=seed;Length=length;Alpha=matchscore;Beta=mismatch;Delta=gap;maxwraplength=wrap;Min_Distance_Window=window;
    Rows=17;Wrapend=Reportmin=Maxrealrow=Maxrow=Maxcol=Maxscore=-73;
    for(int r=0;r<ROWS;r++) {rowptrs[r]=matrix[r];for(int c=0;c<COLS;c++) matrix[r][c]=(int)(rnd()%500)-250;}
    for(int c=0;c<COLS;c++){EC[c]=(unsigned char)(rnd()%256);Up[c]=(int)(rnd()%500)-250;Diag[c]=(int)(rnd()%500)-250;}
    const char *alphabet="ACGTNRY";
    for(int i=0;i<270;i++) Sequence[i]=alphabet[rnd()%7];
    for(int c=0;c<COLS;c++) Consensus.pattern[c]=alphabet[rnd()%7];
    if(style==0) for(int i=1;i<=Length;i++) Sequence[i]='A';
    if(style==1) for(int i=1;i<=Length;i++) Sequence[i]=alphabet[(i-1)%4];
    if(style==2) for(int i=1;i<=Length;i++) Sequence[i]='N';
    if(style==3) for(int i=1;i<=Length;i++) Sequence[i]=Consensus.pattern[(i-1)%size];
    if(style==4) {
        for(int i=1;i<=Length;i++) Sequence[i]=alphabet[(i-1)%size%4];
        for(int i=start+1;i<=Length;i++) if(i%5==0) Sequence[i]='N';
    }
    for(int i=0;i<65536;i++) SM[i]=mismatch;
    for(int i=0;i<4;i++) SM[alphabet[i]*256+alphabet[i]]=matchscore;
}
struct final_state {int matrix[ROWS][COLS],up[COLS],diag[COLS];unsigned char ec[COLS];int globals[6];double rows;};
static struct final_state result;
static struct alignment {int length,score;char textprime[4096],textsecnd[4096];int indexprime[4096],indexsecnd[4096];} AlignPair,saved_alignment;
static double Copynumber,saved_copies;
/* state read only by the narrowband traceback of the same module, which this test does not call */
#define LOCAL 1
#define GLOBAL 0
static int bandcenter_storage[ROWS], *Bandcenter=bandcenter_storage, Maxrealcol, ldong;

static void save_final(struct final_state *out) {
    memcpy(out->matrix,matrix,sizeof matrix);memcpy(out->up,Up,sizeof out->up);memcpy(out->diag,Diag,sizeof out->diag);memcpy(out->ec,EC,sizeof EC);
    out->globals[0]=Wrapend;out->globals[1]=Reportmin;out->globals[2]=Maxrealrow;out->globals[3]=Maxrow;out->globals[4]=Maxcol;out->globals[5]=Maxscore;out->rows=Rows;
}
'''

MAIN = r'''
int main(void) {
    projection=PROJECTION;
    const int scores[][3]={{2,-5,-7},{2,-7,-7},{1,-1,-1},{1,0,0},{0,0,0},{0,-5,-7},{7,-1,0},{5,3,-2},{-2,-3,-1},{2,-5,2},{0,0,1},{100,-100,-1000}};
    for(int size=1;size<=25;size++) for(int cons=0;cons<=1;cons++)
    for(unsigned int scoring=0;scoring<sizeof scores/sizeof scores[0];scoring++)
    for(int style=0;style<7;style++) for(int edge=0;edge<8;edge++) {
        int start=size+(edge==0?0:edge==1?1:edge==2?4:edge==3?20:35);
        int length=edge==4?start:edge==5?start+1:start+43;
        int wrap=edge==6?1:edge==7?2:edge==0?3:70;
        int window=edge%3==0?0:edge%3==1?size:30;
        unsigned int seed=0x184dcf29u+testcase*7919u;
        candidate=0;ntrace=trace_index=0;
        setup(seed,size,style,start,length,wrap,scores[scoring][0],scores[scoring][1],scores[scoring][2],window);
        original_newwrap(start,size,cons);save_final(&result);
        int traceback=scores[scoring][0]>0 && scores[scoring][1]<=0 && scores[scoring][2]<0 && Maxscore>0;
        if(traceback) {traceback_cases++;memset(&AlignPair,0,sizeof AlignPair);get_pair_alignment_with_copynumber(size);saved_alignment=AlignPair;saved_copies=Copynumber;}
        candidate=1;trace_index=0;
        setup(seed,size,style,start,length,wrap,scores[scoring][0],scores[scoring][1],scores[scoring][2],window);
        CANDIDATE_FUNCTION(start,size,cons);
        struct final_state actual;memset(&actual,0,sizeof actual);save_final(&actual);
        int matrix_diff=0;
        if(projection && Maxscore>0 && Maxrow<ROWS) {
            for(int rr=0;rr<=Maxrow;rr++) if(memcmp(result.matrix[rr],actual.matrix[rr],size*sizeof(int))) matrix_diff=1;
        } else matrix_diff=memcmp(result.matrix,actual.matrix,sizeof result.matrix);
        if(traceback) {memset(&AlignPair,0,sizeof AlignPair);get_pair_alignment_with_copynumber(size);if(memcmp(&saved_alignment,&AlignPair,sizeof AlignPair)||saved_copies!=Copynumber) {fprintf(stderr,"traceback differs case %u\n",testcase);return 1;}}
        if(trace_index!=ntrace || matrix_diff || memcmp(result.up,actual.up,sizeof result.up) || memcmp(result.diag,actual.diag,sizeof result.diag) || memcmp(result.ec,actual.ec,sizeof result.ec) || memcmp(result.globals+(projection?1:0),actual.globals+(projection?1:0),(6-(projection?1:0))*sizeof(int)) || (!projection && result.rows!=actual.rows)) {
            fprintf(stderr,"final case %u differs: size=%d consensus=%d scoring=%u style=%d edge=%d\n",testcase,size,cons,scoring,style,edge);return 1;
        }
        testcase++;
    }
    printf("%u %s cases; %u real tracebacks; phase rows %u/%u/%u; AVX2 calls %lu\n",testcase,projection?"endpoint/live-matrix projection":"full-state",traceback_cases,phase_seen[0],phase_seen[1],phase_seen[2],avx2_calls);
    return 0;
}
'''


def kernel_test(directory, original=True, projection=False):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    reference = original_block() if original else fallback_block()
    reference = reference.replace('void newwrap(', 'void original_newwrap(', 1)
    candidate_text = traced(MODULE.read_text())
    # Count calls answered by the AVX2 rows (the similarity matrix follows Alpha/Beta).
    needle = '\t\tReportmin=0;\n\t\t/* A positive accurate maximum defines all reported endpoints. Preserve'
    assert candidate_text.count(needle) == 1
    candidate_text = candidate_text.replace(needle, '\t\tavx2_calls++;\n' + needle)
    source = ('#define TRF_TEST_FULLSTATE 1\n' + SUPPORT + define_guards() + traced(reference) + '\n' +
              candidate_text + '\n#include "alignment/pair_alignment.c"\n' + MAIN)
    source = source.replace('CANDIDATE_FUNCTION', 'newwrap' if projection else 'newwrap_fullstate').replace('PROJECTION', str(int(projection)))
    path = directory / 'kernel.c'
    path.write_text(source)
    binary = directory / 'kernel'
    command = [*cc(), *CFLAGS, '-I', str(SRC), str(path), '-o', str(binary)]
    extra = shlex.split(os.environ.get('TRF_KERNEL_CFLAGS', ''))
    subprocess.run(command + extra, check=True)
    result = subprocess.run([str(binary)], capture_output=True, text=True)
    (directory / 'stdout').write_text(result.stdout)
    (directory / 'stderr').write_text(result.stderr)
    if result.returncode or result.stderr:
        raise RuntimeError(f'kernel failed (exit {result.returncode}); see {directory / "stderr"}')
    print(result.stdout.strip())


class NewwrapTest(unittest.TestCase):
    def test_fallback_is_the_original_function(self):
        self.assertEqual(fallback_block(), original_block())

    def test_complete_state(self):
        with tempfile.TemporaryDirectory(prefix='trf-newwrap-kernel-') as directory:
            kernel_test(directory)

    def test_traceback_projection(self):
        with tempfile.TemporaryDirectory(prefix='trf-newwrap-projection-') as directory:
            kernel_test(directory, projection=True)

    def test_isolated_end_to_end(self):
        with tempfile.TemporaryDirectory(prefix='trf-newwrap-exact-') as directory:
            directory = Path(directory)
            binary = build(directory / 'build')
            for case in exact_cases():
                if case.group == 'long':
                    continue
                with self.subTest(case=case.name):
                    left = directory / case.name / 'reference'
                    right = directory / case.name / 'candidate'
                    reference = execute(oracle_dir() / 'trf', case, left, 120)
                    candidate = execute(binary, case, right, 120)
                    self.assertEqual(reference['exit'], case.expected_exit)
                    self.assertEqual(candidate['exit'], case.expected_exit)
                    self.assertEqual(compare(reference, candidate, left, right), [])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build-dir', type=Path)
    parser.add_argument('--kernel-dir', type=Path)
    parser.add_argument('--kernel-only', action='store_true')
    parser.add_argument('--oracle', type=Path,
                        help='directory containing source.tar and trf of the original; overrides TRF_ORACLE_DIR')
    options, remaining = parser.parse_known_args()
    if options.oracle:
        os.environ['TRF_ORACLE_DIR'] = str(options.oracle.resolve())
    if options.kernel_only:
        with tempfile.TemporaryDirectory(prefix='trf-newwrap-kernel-') as directory:
            directory = options.kernel_dir.resolve() if options.kernel_dir else Path(directory)
            kernel_test(directory, original=False)
            kernel_test(directory / 'projection', original=False, projection=True)
        raise SystemExit(0)
    require_oracle()
    if options.build_dir:
        print(build(options.build_dir.resolve()))
    elif options.kernel_dir:
        kernel_test(options.kernel_dir.resolve())
        kernel_test(options.kernel_dir.resolve() / 'projection', projection=True)
    else:
        unittest.main(argv=[__file__, *remaining])
