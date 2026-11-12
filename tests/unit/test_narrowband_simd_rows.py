#!/usr/bin/env python3
"""Compare the SIMD row helpers with the ordered scalar finite-sentinel recurrence."""
import argparse
import os
from pathlib import Path
import platform
import subprocess
import sys
import tempfile

sys.path[:0] = [str(Path(__file__).resolve().parents[1] / name) for name in ('support', 'integration')]
from trf_test_support import SRC, cc, require_sanitizers
from test_narrowband_simd import preamble_with_beta


def main():
    with tempfile.TemporaryDirectory(prefix='trf-narrowband-rows-') as temporary:
        run(temporary)


def run(temporary):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path,
                        help='new directory for the harness and its output; temporary if omitted')
    args = parser.parse_args()
    if platform.machine() not in ("x86_64", "i386", "i686"):
        print("SKIP: SSE4.1 row helper is x86 specific")
        return
    # The harness supplies the program state and compiles the complete vector
    # kernel module; the row helpers under test are static functions of it.
    # The AVX2 run loops own their row bookkeeping and are covered by the
    # complete-state comparison in test_narrowband_simd.py.
    program = preamble_with_beta() + r'''
#include <limits.h>
#include <smmintrin.h>
#include "alignment/narrowband_simd.c"
#undef max2
#undef max3
#define max2(a,b) ((a)>=(b)?(a):(b))
#define max3(a,b,c) max2(max2(a,b),c)
__attribute__((target("sse4.1"),noinline))
static void invoke_ramp(int *current,const int *previous,const unsigned char *pattern,
                   int size,int c,int w,int direction,int gate,int character,
                   struct trf_narrowband_row *out) {
    if(direction>0) trf_narrowband_scan_forward(current,previous,pattern,size,c,w,gate,character,out);
    else trf_narrowband_scan_backward(current,previous,pattern,size,c,w,gate,character,out);
}
__attribute__((target("sse4.1"),noinline))
static void invoke(int *current,const int *previous,const unsigned char *pattern,
                   int size,int c,int w,int direction,int gate,int character,
                   struct trf_narrowband_row *out) {
    trf_narrowband_scan(current,previous,pattern,size,c,w,direction,gate,character,out);
}
static uint32_t seed=0xb793a1d5;
static uint32_t next(void) {seed=1664525u*seed+1013904223u;return seed;}
static int canonical(int c) {return c=='A'||c=='C'||c=='G'||c=='T';}
static void baseline(int *current,const int *previous,const unsigned char *pattern,
                     int size,int c,int w,int direction,int gate,int character,
                     struct trf_narrowband_row *out) {
    int left=-1000,rowbest=-1,score=-1000,last=0,matchcol=-71;
    for(int t=0;t<=2*w;t++) {
        int i=direction>0?t:2*w-t;
        int sim=(canonical(character)&&pattern[c]==character)?Alpha:Beta;
        int diag=previous[i]+sim;
        int value=max3(0,diag,left);
        if(t<2*w)value=max2(value,previous[i+direction]+Delta);
        left=value+Delta;
        if(gate&&value==0)value=left=-1000;
        current[i]=value;
        if(value>=score){score=value;last=t;}
        if(value>rowbest){rowbest=value;matchcol=(value==diag&&sim==Alpha)?c:-2;}
        c=(c+direction+size)%size;
    }
    out->score=score;out->last=last;out->matchcol=matchcol;out->dead=gate&&score==-1000;
}
int main(void) {
    maxwraplength=100;
    if(!__builtin_cpu_supports("sse4.1")){puts("SKIP: SSE4.1 unavailable");return 0;}
    const int alphas[]={0,1,2,999,1000,1001},betas[]={0,-1,-7,-999,-1000,-1001};
    const int gaps[]={0,-1,-7,-999,-1000,-1001},widths[]={1,2,3,4,5,6,7,75};
    const int values[]={-1000,-999,-1,0,1,7,999,1000,1001};
    const unsigned char chars[]="ANT";
    int previous[151],actual[151],expected[151];unsigned char pattern[151];
    unsigned int cases=0,dead=0,zero=0,ramp=0;
    for(unsigned a=0;a<6;a++)for(unsigned b=0;b<6;b++)for(unsigned d=0;d<6;d++)
    for(unsigned wi=0;wi<8;wi++)for(int direction=-1;direction<=1;direction+=2)
    for(int gate=0;gate<=1;gate++)for(int style=0;style<3;style++)for(int ch=0;ch<3;ch++) {
        int w=widths[wi],size=2*w+1,c=style==0?0:style==1?size-1:size/2;
        Alpha=alphas[a];Beta=betas[b];Delta=gaps[d];
        /* cells beyond the band hold arbitrary values and must be preserved */
        for(int i=0;i<151;i++) previous[i]=values[next()%9];
        for(int i=0;i<size;i++) {
            previous[i]=style==0?0:style==1?-1000:values[next()%9];
            pattern[i]="ACGTNRY"[next()%7];
        }
        struct trf_narrowband_row wanted,got;
        for(int i=0;i<151;i++) expected[i]=0x5a5a0000+i;
        baseline(expected,previous,pattern,size,c,w,direction,gate,chars[ch],&wanted);
        for(int path=0;path<3;path++) {
        for(int i=0;i<151;i++) actual[i]=0x5a5a0000+i;
        if(path==2) {
            /* isolated rows: previous cells are at most 1001 here */
            if(!trf_narrowband_ramp_domain(w,w) || 1001>trf_narrowband_packed_limit()) continue;
            invoke_ramp(actual,previous,pattern,size,c,w,direction,gate,chars[ch],&got);
            ramp++;
        }
        else if(path) trf_narrowband_scan_unpacked(actual,previous,pattern,size,c,w,direction,gate,chars[ch],&got);
        else invoke(actual,previous,pattern,size,c,w,direction,gate,chars[ch],&got);
        if(memcmp(expected,actual,sizeof actual) || wanted.score!=got.score || wanted.last!=got.last || wanted.dead!=got.dead || (!wanted.dead&&wanted.matchcol!=got.matchcol)) {
            fprintf(stderr,"path=%d ",path);
            fprintf(stderr,"case=%u A=%d B=%d D=%d width=%d dir=%d gate=%d style=%d char=%c expected=%d/%d/%d/%d actual=%d/%d/%d/%d\n",cases,Alpha,Beta,Delta,w,direction,gate,style,chars[ch],wanted.score,wanted.last,wanted.matchcol,wanted.dead,got.score,got.last,got.matchcol,got.dead);return 1;
        }
        cases++;dead+=got.dead;zero+=got.score==0;
        }
    }
    const int packed_edges[]={8388607,8388608,INT_MAX};
    maxwraplength=1;Beta=Delta=0;
    for(int edge=0;edge<3;edge++)for(int w=1;w<=75;w++)
    for(int direction=-1;direction<=1;direction+=2)for(int gate=0;gate<=1;gate++) {
        Alpha=packed_edges[edge];int size=2*w+1,c=direction>0?size-1:0;
        for(int i=0;i<size;i++){previous[i]=0;pattern[i]='A';}
        struct trf_narrowband_row wanted,got;
        baseline(expected,previous,pattern,size,c,w,direction,gate,'A',&wanted);
        invoke(actual,previous,pattern,size,c,w,direction,gate,'A',&got);
        if(memcmp(expected,actual,size*sizeof(int)) || wanted.score!=got.score || wanted.last!=got.last || wanted.matchcol!=got.matchcol || wanted.dead!=got.dead) {
            fprintf(stderr,"packed threshold mismatch A=%d width=%d dir=%d gate=%d\n",Alpha,w,direction,gate);return 2;
        }
        cases++;
    }
    printf("%u SIMD ordered-row cases matched (%u through the ramp kernels); %u all-killed; %u surviving-zero maxima\n",cases,ramp,dead,zero);
    if(!ramp) return 3;
    return 0;
}
'''
    require_sanitizers('-fsanitize=address,undefined')
    output = (args.output or Path(temporary) / 'work').resolve()
    output.mkdir(parents=True, exist_ok=False)
    (output/'rows.c').write_text(program)
    subprocess.run([*cc(),'-O1','-g','-fsanitize=address,undefined','-fno-sanitize-recover=all','-I',str(SRC),str(output/'rows.c'),'-o',str(output/'rows')],check=True)
    result = subprocess.run([str(output/'rows')],capture_output=True,text=True,env={**os.environ,'ASAN_OPTIONS':'detect_leaks=0'})
    (output/'stdout').write_text(result.stdout)
    (output/'stderr').write_text(result.stderr)
    print(result.stdout,end='')
    if result.stderr:print(result.stderr,end='')
    result.check_returncode()


if __name__ == '__main__':
    main()
