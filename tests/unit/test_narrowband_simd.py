#!/usr/bin/env python3
"""Compare all global state of the vector narrowband kernels with the original function.

The similarity matrix follows init_sm: only canonical self-pairs match.
The original function and the modules src/alignment/narrowband_scalar.c,
narrowband_simd.c, and narrowband.c are compiled in one translation unit
against the harness's copy of the program state, with no global ISA flag.
All scratch and padding is poisoned. Counters confirm that the vector rows,
the scalar kernel, and the unencoded scan were all exercised.
"""
import argparse
import os
from pathlib import Path
import subprocess
import sys
import tempfile

sys.path[:0] = [str(Path(__file__).resolve().parents[1] / name) for name in ('support', 'integration')]
from trf_test_support import SRC, cc, original_source, require_oracle
from test_narrowband_state import SCALAR, DRIVER, block_bounds, instrument, preamble

SIMD=SRC/'alignment/narrowband_simd.c'
DISPATCH=SRC/'alignment/narrowband.c'
VECTOR='__attribute__((target("sse4.1")))\nvoid trf_narrowband_vector('


def _shift_branches(text):
    """The two scalar band-shift branches (backward, forward), whitespace-normalised."""
    import re
    out=[]
    for start in (r'/\* band shifts left \*/', r'if\(k>=1\)\s+/\* band shifts right \*/'):
        begin=re.search(start,text)
        assert begin, start
        end=text.index('lastmatchatmax_col+size',begin.end())
        end=text.rindex('\n',0,end)
        out.append(' '.join(text[begin.end():end].split()))
    return out


def check_scalar_clone():
    """The scalar band-shift branches kept inside the vector function must stay
    identical to the scalar kernel. (The surrounding row bookkeeping
    of the vector function intentionally differs: band-center prediction and
    selections instead of branches; the complete-state comparison covers it.)"""
    scalar=SCALAR.read_text()
    candidate=SIMD.read_text()
    clone=candidate[candidate.index(VECTOR):]
    reference=_shift_branches(scalar[scalar.index('void trf_narrowband_scalar('):])
    copied=_shift_branches(clone)
    assert copied==reference, 'scalar rows of the vector function differ from the scalar kernel'


def preamble_with_beta():
    text=preamble()
    for old,new in (('int alpha, delta,','int alpha, beta, delta,'),
                    ('#define Delta state.delta','#define Beta state.beta\n#define Delta state.delta')):
        assert text.count(old)==1
        text=text.replace(old,new)
    return text


def run(directory, sanitize=False, force_scalar=False, trials=12000, portable=False, variant='avx2'):
    check_scalar_clone()
    directory=Path(directory)
    directory.mkdir(parents=True,exist_ok=True)
    original=original_source('src/tr30dat.c')
    start,end=block_bounds(original)
    original=instrument(original[start:end]).replace('void narrowbandwrap(', 'void original_narrowbandwrap(')
    undef='\n'.join('#undef '+name for name in (
        'test_trace_and_backwards_maxscore','test_maxrowscore_with_match',
        'test_maxrowscore_without_match','test_trace_and_forward_maxscore'))
    # Instrument each function separately; shape counters participate in the
    # complete state comparison, whereas dispatch counters are observation only.
    candidate=SIMD.read_text()
    vector_at=candidate.index(VECTOR)
    candidate=instrument(SCALAR.read_text())+candidate[:vector_at]+instrument(candidate[vector_at:])+DISPATCH.read_text()
    candidate=candidate.replace('    __m128i zero=_mm_setzero_si128()',
                                '    ++vector_rows[direction>0];\n    __m128i zero=_mm_setzero_si128()')
    # AVX2 run loops: same observation counters as the instrumented row loop.
    for phase,sign in ((1,2),(0,0)):
        name='trf_narrowband_run_forward_avx2' if phase else 'trf_narrowband_run_backward_avx2'
        at=candidate.index('static void '+name+'(')
        needle='        rows++;\n'
        pos=candidate.index(needle,at)
        candidate=(candidate[:pos]+needle+
                   f'        state.coverage[{phase}][{sign}]++; state.coverage[{phase}][3]+=recenter; ++vector_rows[{phase}]; ++avx2_rows_seen;\n'+
                   candidate[pos+len(needle):])
    candidate=candidate.replace('        trf_narrowband_scan_unpacked(current,previous,', '        ++unpacked_rows;\n        trf_narrowband_scan_unpacked(current,previous,')
    candidate=candidate.replace('    trf_narrowband_scalar(start,size,',
                                '    ++scalar_calls;\n    trf_narrowband_scalar(start,size,')
    preamble=preamble_with_beta()
    preamble+='\nstatic unsigned long vector_rows[2],scalar_calls,unpacked_rows,avx2_rows_seen;\n'
    if variant=='sse-ramp': preamble+='\n#define TRF_NARROWBAND_DISABLE_AVX2 1\n'
    if variant=='packed': preamble+='\n#define TRF_NARROWBAND_DISABLE_AVX2 1\n#define TRF_NARROWBAND_DISABLE_RAMP 1\n'
    if force_scalar: preamble+='\n#define TRF_NARROWBAND_FORCE_SCALAR 1\n'
    if portable: preamble+='\n#define TRF_NARROWBAND_DISABLE_SIMD 1\n'
    driver=DRIVER.replace('trial<6000',f'trial<{trials}')
    driver=driver.replace('Alpha=1+(trial%4);','Alpha=(trial%97==0)?100000:(trial%29==0)?0:(trial%31==0?1001:1+(trial%4)); Beta=(trial%41==0)?0:(trial%43==0?-1001:-1-(trial%7));')
    driver=driver.replace('Delta=(trial%17==0)?0:', 'Delta=(trial%47==0)?-1001:(trial%17==0)?0:')
    driver=driver.replace('Rows=trial%113;', 'if(Alpha==0) Reportmin=100000;\n        Rows=trial%113;')
    driver=driver.replace('?Alpha:-1-(trial%7);','?Alpha:Beta;')
    driver=driver.replace('puts("6000 complete-state original/candidate cases matched");',
        f'printf("{trials} complete-state original/candidate cases matched ({variant}); vector rows=%lu/%lu (avx2 run rows=%lu) scalar calls=%lu unencoded rows=%lu\\n",vector_rows[0],vector_rows[1],avx2_rows_seen,scalar_calls,unpacked_rows);\n'
        +('if(vector_rows[0] || vector_rows[1] || !scalar_calls) return 3;' if force_scalar or portable else 'if(!vector_rows[0] || !vector_rows[1] || !scalar_calls || !unpacked_rows) return 3;'))
    source=directory/'state.c'
    source.write_text(preamble+original+'\n'+undef+'\n'+candidate+driver)
    binary=directory/'state'
    flags=['-O1','-g','-fsanitize=address,undefined','-fno-sanitize-recover=all'] if sanitize else ['-O2','-g']
    subprocess.run([*cc(),*flags,'-I',str(SRC),str(source),'-o',str(binary)],check=True)
    subprocess.run([str(binary)],check=True,env=dict(os.environ,ASAN_OPTIONS='detect_leaks=0'))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sanitize',action='store_true')
    parser.add_argument('--force-scalar',action='store_true')
    parser.add_argument('--portable',action='store_true')
    parser.add_argument('--trials',type=int,default=12000)
    parser.add_argument('--build-dir',type=Path)
    parser.add_argument('--variant',choices=('avx2','sse-ramp','packed'),default='avx2',
                        help='row kernels: AVX2 run loops, SSE4.1 ramp rows, or the packed SSE4.1 scan')
    args=parser.parse_args()
    require_oracle('source.tar')
    if args.build_dir: run(args.build_dir,args.sanitize,args.force_scalar,args.trials,args.portable,args.variant)
    else:
        with tempfile.TemporaryDirectory(prefix='trf-narrowband-simd-') as directory:
            run(directory,args.sanitize,args.force_scalar,args.trials,args.portable,args.variant)
