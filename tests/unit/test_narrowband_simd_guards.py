#!/usr/bin/env python3
"""Probe the boundaries of the vector-kernel dispatch test under UBSan.

The harness links src/alignment/narrowband_simd.c and calls
trf_narrowband_vector_domain with the program's own scoring globals; no
alignment matrix is allocated."""
import argparse
from pathlib import Path
import platform
import subprocess
import sys
import tempfile

sys.path[:0] = [str(Path(__file__).resolve().parents[1] / name) for name in ('support', 'integration')]
from trf_test_support import SRC, cc, require_sanitizers

MODULES = ('alignment/narrowband_simd.c', 'core/globals.c', 'core/report.c')


def main():
    with tempfile.TemporaryDirectory(prefix='trf-narrowband-guards-') as temporary:
        run(temporary)


def run(temporary):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path,
                        help='new directory for the harness and its output; temporary if omitted')
    args = parser.parse_args()
    if platform.machine() not in ('x86_64', 'i386', 'i686'):
        print('SKIP: dispatch helper is x86 specific')
        return
    require_sanitizers('-fsanitize=undefined')
    output = (args.output or Path(temporary) / 'work').resolve()
    output.mkdir(parents=True, exist_ok=False)
    program = r'''
#include <assert.h>
#include <stdint.h>
#include <limits.h>
#include <stdio.h>
#include "core/globals.h"
#include "alignment/narrowband_simd.h"
/* widest band radius of a row of MAXBANDWIDTH cells, and the widest band */
#define RADIUS ((MAXBANDWIDTH-1)/2)
#define WIDEST (2*RADIUS+1)
static int available;
static void defaults(void) { Alpha=2;Beta=-5;Delta=-7;maxwraplength=100; }
static void yes(int size,int wb,int wf) {
    assert(trf_narrowband_vector_domain(size,size,wb,wf)==available);
}
static void no(int start,int size,int wb,int wf) {
    assert(!trf_narrowband_vector_domain(start,size,wb,wf));
}
int main(void) {
    available=__builtin_cpu_supports("sse4.1")!=0;
    defaults(); yes(3,1,1);yes(WIDEST,RADIUS,RADIUS);yes(4000,RADIUS,RADIUS);
    no(4001,4001,1,1);no(2,2,1,1);no(2,3,1,1);
    no(100,100,0,1);no(100,100,1,0);no(100,100,-1,1);
    no(100,100,INT_MAX,1);no(100,100,1,INT_MAX);
    no(4000,4000,RADIUS+1,1);no(4000,4000,1,RADIUS+1);no(100,100,50,1);
    maxwraplength=0;no(3,3,1,1);
    defaults();Alpha=-1;no(3,3,1,1);
    defaults();Beta=1;no(3,3,1,1);
    defaults();Delta=1;no(3,3,1,1);
    defaults();Delta=INT_MIN;no(3,3,1,1);
    defaults();Beta=INT_MIN;no(3,3,1,1);
    defaults();Alpha=INT_MAX;maxwraplength=1;yes(3,1,1);
    maxwraplength=2;no(3,3,1,1);
    Alpha=INT_MAX/2;yes(3,1,1);
    Alpha=INT_MAX/2+1;no(3,3,1,1);
    Alpha=0;maxwraplength=UINT_MAX;yes(3,1,1);
    defaults();Delta=-(INT_MAX/8);yes(3,1,1);
    Delta=-(INT_MAX/8+1);no(3,3,1,1);
    defaults();Delta=0;Beta=-(INT_MAX-999);yes(3,1,1);
    Beta=-(INT_MAX-998);no(3,3,1,1);
    /* (wb+wf)*d+d may reach, but not exceed, 2^31 */
    defaults();Delta=-(int)((-(int64_t)INT_MIN)/WIDEST);yes(WIDEST,RADIUS,RADIUS);
    --Delta;no(WIDEST,WIDEST,RADIUS,RADIUS);
    puts("SIMD dispatch boundary assertions passed under UBSan");
    return 0;
}
'''
    (output/'guards.c').write_text(program)
    subprocess.run([*cc(),'-O1','-g','-fsanitize=undefined','-fno-sanitize-recover=all','-I',str(SRC),
                    str(output/'guards.c'),*(str(SRC/name) for name in MODULES),'-lm',
                    '-o',str(output/'guards')],check=True)
    result = subprocess.run([str(output/'guards')],capture_output=True,text=True)
    (output/'stdout').write_text(result.stdout)
    (output/'stderr').write_text(result.stderr)
    print(result.stdout,end='')
    if result.stderr:
        print(result.stderr,end='')
    result.check_returncode()


if __name__ == '__main__':
    main()
