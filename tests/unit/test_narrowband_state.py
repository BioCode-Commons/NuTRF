#!/usr/bin/env python3
"""Compare all narrowband global state with the original function.

The harness holds the program state in one structure and compiles two functions
against it: narrowbandwrap from the original TRF sources, and the scalar kernel
module src/alignment/narrowband_scalar.c with its includes of the shared
headers disabled. The executable runs both functions from identical seeded state,
including poisoned unused rows/cells/scratch, positive/zero/negative gaps,
zero-iteration phases, and both larger and smaller forward bands. It compares every byte of the
matrix, Bandcenter, Diag, Up, EC, input/scoring data, and endpoint/counter state.
Use --sanitize for ASan/UBSan on this isolated kernel (excluding other kernels).
"""
import argparse
import os
from pathlib import Path
import subprocess
import sys
import tempfile

sys.path[:0] = [str(Path(__file__).resolve().parents[1] / name) for name in ('support', 'integration')]
from trf_test_support import SRC, cc, define_guards, original_source, require_oracle

SCALAR = SRC / 'alignment/narrowband_scalar.c'
# Headers that belong to the narrowband kernels stay in effect; the harness
# replaces everything the other headers declare.
KERNEL_HEADERS = ('alignment/narrowband.h', 'alignment/narrowband_scalar.h',
                  'alignment/narrowband_simd.h')


def block_bounds(source):
    """The four narrowband macros and complete function; no neighboring macros."""
    start = source.index('/***************************** narrowbandwrap()  macros')
    start = source.rfind('/*******************************************************************/', 0, start)
    following = source.index('/***************************** globalnarrowbandwrap()  macros')
    end = source.rfind('}', start, following) + 1
    return start, end


PREAMBLE = r'''
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdarg.h>
#define MAXBANDWIDTH 151
#define MAXR 128
#define MAXSEQ 4096
#define WITHCONSENSUS 1
#define FALSE 0
#define TRUE 1
#define max2(a,b) ((a>=b)?a:b)
#define max3(a,b,c) ((a>=b)?((a>=c)?a:c):((b>=c)?b:c))
#define max4(a,b,c,d) ((a>=b)?((a>=c)?((a>=d)?a:d):((c>=d)?c:d)):((b>=c)?((b>=d)?b:d):((c>=d)?c:d)))
#define match(a,b) SM[256*(a)+(b)]
struct state {
    int scores[MAXR+1][MAXBANDWIDTH+1];
    int centers[MAXR+1], diag[MAXBANDWIDTH+1], up[MAXBANDWIDTH+1];
    unsigned char ec[4001], sequence[MAXSEQ+1];
    struct { unsigned char pattern[4001]; } consensus;
    int sm[65536];
    int alpha, delta, length, ldong, reportmin, criteria_print;
    unsigned int maxwrap;
    double rows;
    int wrapend, maxrealrow, maxrow, maxcol, maxscore;
    unsigned int coverage[2][4];
};
static struct state state, before, expected;
static int *matrix_rows[MAXR+1];
#define S matrix_rows
#define Bandcenter state.centers
#define Diag state.diag
#define Up state.up
#define EC state.ec
#define Sequence state.sequence
#define Consensus state.consensus
#define SM state.sm
#define Alpha state.alpha
#define Delta state.delta
#define Length state.length
#define ldong state.ldong
#define Reportmin state.reportmin
#define Criteria_print state.criteria_print
#define maxwraplength state.maxwrap
#define Rows state.rows
#define Wrapend state.wrapend
#define Maxrealrow state.maxrealrow
#define Maxrow state.maxrow
#define Maxcol state.maxcol
#define Maxscore state.maxscore
static void trf_message(const char *format, ...) { (void)format; }
'''


def preamble():
    return PREAMBLE + define_guards(KERNEL_HEADERS)

DRIVER = r'''
static uint32_t seed=0x5a727646;
static uint32_t next(void) { seed=1664525u*seed+1013904223u; return seed>>8; }
static const unsigned char alphabet[]="ACGTNNRY";
int main(void) {
    unsigned long coverage[2][4]={{0}};
    for (int row=0; row<=MAXR; ++row) matrix_rows[row]=state.scores[row];
    for (int trial=0; trial<6000; ++trial) {
        /* Poison all untouched state and padding, then constrain live inputs. */
        memset(&state,0x35+(trial%19),sizeof state);
        int size=21+(next()%380);
        if (trial%11==0) size=1999+(trial%3);
        if (trial%19==0) size=3+(trial%18);
        int maxwidth=(size-1)/2;
        if (maxwidth>75) maxwidth=75;
        int backward=1+(next()%maxwidth), forward=1+(next()%maxwidth);
        if (trial%7==0) forward=backward;
        int start=size+1+(next()%160);
        int option=trial%2, tuplesize=1+trial%7;
        maxwraplength=1+(next()%MAXR);
        Length=start+(next()%180);
        if (trial%17==0) Length=start;
        ldong=(trial%5==0);
        Alpha=1+(trial%4);
        Delta=(trial%17==0)?0:(trial%23==0?1:-1-(trial%9));
        Reportmin=(trial%3==0)?100000:(trial%2?1:50);
        Criteria_print=trial%2;
        if (trial%37==0) { maxwraplength=0; ldong=1; }
        if (trial%31==0) Length=0; /* No forward rows: preserve initialization. */
        /* Zero reporting threshold is defined here: homopolymer calls have
         * at least one forward row and a positive final score. */
        if (trial%4==1 && Length>start && maxwraplength>0) Reportmin=0;
        Rows=trial%113;
        memset(state.coverage,0,sizeof state.coverage);
        for (int a=0;a<256;++a) for (int b=0;b<256;++b)
            SM[256*a+b]=(a==b && (a=='A'||a=='C'||a=='G'||a=='T'))?Alpha:-1-(trial%7);
        for (int i=0;i<=MAXSEQ;++i) Sequence[i]=alphabet[next()%8];
        for (int i=0;i<size;++i) {
            Consensus.pattern[i]=alphabet[next()%8];
            /* Periodic, random, shifted, and substitution/indel-like streams. */
            if (trial%4!=0) {
                unsigned char c=alphabet[next()%4];
                if (trial%4==1) c='A';
                Consensus.pattern[i]=c;
                for (int j=i+1;j<=Length;j+=size) Sequence[j]=c;
            }
        }
        if (trial%4==3) {
            for (int j=start-size;j<=Length;++j)
                if (j>0 && j%19==0) Sequence[j]=alphabet[next()%8];
        }
        before=state;
        original_narrowbandwrap(start,size,backward,forward,option,tuplesize);
        expected=state;
        state=before;
        narrowbandwrap(start,size,backward,forward,option,tuplesize);
        if (memcmp(&state,&expected,sizeof state)) {
            unsigned char *a=(unsigned char *)&state, *b=(unsigned char *)&expected;
            size_t at=0; while(a[at]==b[at]) ++at;
            fprintf(stderr,"trial=%d size=%d widths=%d/%d option=%d tuple=%d ldong=%d byte=%zu actual=%u expected=%u\n",
                    trial,size,backward,forward,option,tuplesize,ldong,at,a[at],b[at]);
            fprintf(stderr,"maxwrap=%u rows=%g/%g\n",maxwraplength,Rows,expected.rows);
            return 1;
        }
        for(int phase=0;phase<2;++phase) for(int i=0;i<4;++i)
            coverage[phase][i]+=state.coverage[phase][i];
    }
    for(int phase=0;phase<2;++phase) {
        printf("phase=%d k-negative=%lu k-zero=%lu k-positive=%lu recentered=%lu\n",
               phase,coverage[phase][0],coverage[phase][1],coverage[phase][2],coverage[phase][3]);
        for(int i=0;i<4;++i) if(!coverage[phase][i]) return 2;
    }
    puts("6000 complete-state original/candidate cases matched");
    return 0;
}
'''


def instrument(block):
    # Observation-only counters force evidence for both k signs, zero, recenter.
    needle = 'if(size-k<=k) k=-(size-k);'
    assert block.count(needle) == 2
    for phase in range(2):
        pos = block.index(needle)
        # The marker prevents matching this already-instrumented occurrence.
        replacement = ('if(size-k<=k) k= -(size-k);' +
                       f'\nstate.coverage[{phase}][k<0?0:k==0?1:2]++;' +
                       f'\nstate.coverage[{phase}][3]+=(matches_in_diagonal>=tuplesize);')
        block = block[:pos] + block[pos:].replace(needle, replacement, 1)
    return block


def run(directory, sanitize=False, module=SCALAR):
    directory=Path(directory)
    directory.mkdir(parents=True,exist_ok=True)
    original=original_source('src/tr30dat.c')
    start,end=block_bounds(original)
    original=instrument(original[start:end]).replace('void narrowbandwrap(', 'void original_narrowbandwrap(')
    undef='\n'.join('#undef '+name for name in (
        'test_trace_and_backwards_maxscore', 'test_maxrowscore_with_match',
        'test_maxrowscore_without_match', 'test_trace_and_forward_maxscore'))
    source=directory/'state.c'
    driver=DRIVER.replace('        narrowbandwrap(start,size,','        trf_narrowband_scalar(start,size,')
    assert driver!=DRIVER
    source.write_text(preamble()+original+'\n'+undef+'\n'+instrument(Path(module).read_text())+driver)
    binary=directory/'state'
    flags=['-O1','-g','-fsanitize=address,undefined','-fno-sanitize-recover=all'] if sanitize else ['-O2','-g']
    subprocess.run([*cc(),*flags,'-I',str(SRC),str(source),'-o',str(binary)],check=True)
    environment=dict(os.environ,ASAN_OPTIONS='detect_leaks=0')
    subprocess.run([str(binary)],check=True,env=environment)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sanitize',action='store_true')
    parser.add_argument('--build-dir',type=Path)
    parser.add_argument('--module',type=Path,default=SCALAR,help='scalar kernel source file to test')
    args=parser.parse_args()
    require_oracle('source.tar')
    if args.build_dir:
        run(args.build_dir,args.sanitize,args.module)
    else:
        with tempfile.TemporaryDirectory(prefix='trf-narrowband-state-') as directory:
            run(directory,args.sanitize,args.module)
