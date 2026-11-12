#!/usr/bin/env python3
"""Compare original TRF and candidate status, streams, filenames and raw bytes."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import signal
import subprocess
import sys
import threading
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'support'))
from trf_corpus import ROOT, SEED, exact_cases
from trf_test_support import DEFAULT_CANDIDATE, SKIP_EXIT, oracle_dir


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def snapshot(work):
    return {str(path.relative_to(work)): {'size': path.stat().st_size, 'sha256': sha256(path)}
            for path in sorted(work.rglob('*')) if path.is_file() and path.relative_to(work) != Path('input.fasta')}


def execute(binary, case, work, timeout):
    """wait4 measures this child only; streams are retained, never normalized."""
    work.mkdir(parents=True)
    if case.source:
        shutil.copyfile(case.source, work / 'input.fasta')
    else:
        (work / 'input.fasta').write_bytes(case.data)
    arguments = ['trf', *case.arguments()]
    result = {'argv': arguments, 'timeout': False, 'launch_error': None}
    timer = None
    proc = None
    start = time.perf_counter()
    with (work / 'stdout').open('wb') as out, (work / 'stderr').open('wb') as err, \
            (work / 'input.fasta').open('rb') as source:
        try:
            proc = subprocess.Popen(arguments, executable=str(binary), cwd=work,
                                    stdin=source if case.stdin else subprocess.DEVNULL,
                                    stdout=out, stderr=err, start_new_session=True,
                                    env={**os.environ, 'LC_ALL': 'C', 'TZ': 'UTC'})
            lock = threading.Lock()

            def terminate():
                with lock:
                    if proc.returncode is None:
                        try:
                            os.killpg(proc.pid, signal.SIGKILL)
                            result['timeout'] = True
                        except ProcessLookupError:
                            pass

            timer = threading.Timer(timeout, terminate)
            timer.daemon = True
            timer.start()
            _, status, usage = os.wait4(proc.pid, 0)
            with lock:
                proc.returncode = os.waitstatus_to_exitcode(status)
            result.update(exit=proc.returncode, cpu_user_seconds=usage.ru_utime,
                          cpu_system_seconds=usage.ru_stime, peak_rss_kib=usage.ru_maxrss)
        except OSError as error:
            result.update(exit=None, launch_error=str(error))
        finally:
            if timer:
                timer.cancel()
    result['wall_seconds'] = time.perf_counter() - start
    result['files'] = snapshot(work)
    return result


def first_difference(left, right):
    offset = 0
    with left.open('rb') as a, right.open('rb') as b:
        while True:
            x, y = a.read(65536), b.read(65536)
            if x != y:
                n = next((i for i, pair in enumerate(zip(x, y)) if pair[0] != pair[1]), min(len(x), len(y)))
                return {'byte_offset': offset + n, 'reference_hex': x[n:n + 24].hex(),
                        'candidate_hex': y[n:n + 24].hex()}
            if not x:
                return None
            offset += len(x)


def compare(reference, candidate, reference_work, candidate_work):
    differences = []
    for key in ('exit', 'timeout', 'launch_error'):
        if reference[key] != candidate[key]:
            differences.append({'field': key, 'reference': reference[key], 'candidate': candidate[key]})
    left, right = set(reference['files']), set(candidate['files'])
    if left != right:
        differences.append({'field': 'filenames', 'only_reference': sorted(left - right),
                            'only_candidate': sorted(right - left)})
    # Compare actual bytes even when the stored hashes match.
    for name in sorted(left & right):
        difference = first_difference(reference_work / name, candidate_work / name)
        if difference is not None:
            differences.append({'field': 'file_bytes', 'file': name, **difference})
    for label, run in [('reference', reference), ('candidate', candidate)]:
        if run['exit'] is not None and run['exit'] < 0:
            differences.append({'field': 'signal_exit', 'binary': label, 'exit': run['exit']})
    # Both timing out or failing to launch is never a passing differential test.
    if reference['timeout'] or candidate['timeout']:
        differences.append({'field': 'incomplete', 'reason': 'timeout'})
    if reference['launch_error'] or candidate['launch_error']:
        differences.append({'field': 'incomplete', 'reason': 'launch error'})
    return differences


def provenance(reference, candidate):
    manifest = oracle_dir() / 'manifest.json'
    return {'reference': str(reference), 'candidate': str(candidate),
            'reference_sha256': sha256(reference), 'candidate_sha256': sha256(candidate),
            'oracle_manifest': json.loads(manifest.read_text()) if manifest.exists() else None,
            'seed': SEED, 'platform': platform.platform(), 'python': platform.python_version(),
            'cpu_affinity': sorted(os.sched_getaffinity(0)), 'logical_cpu_count': os.cpu_count(),
            'inherited_sanitizer_findings': [
                'Original source: tr30dat.c:690 index -1 in int[151]',
                'Original source: tr30dat.c:2467 signed overflow 46341 * 46341',
                'LeakSanitizer cannot run under this environment ptrace; no leak-clean claim'],
            'environment': {'LC_ALL': 'C', 'TZ': 'UTC'},
            'comparison': 'Identical argv[0]=trf; exit status, stdout, stderr, generated filenames and unmodified file bytes'}


def prepare_output(args):
    output = args.output.resolve()
    if args.overwrite and output.is_dir():
        shutil.rmtree(output)
    output.mkdir(parents=True, exist_ok=False)
    return output


def common_arguments(parser, output_name):
    parser.add_argument('--reference', type=Path, default=oracle_dir() / 'trf',
                        help='original TRF executable (tools/oracle/build_original_trf.sh)')
    parser.add_argument('--candidate', type=Path, default=DEFAULT_CANDIDATE)
    parser.add_argument('--output', type=Path, default=ROOT / 'build' / output_name)
    parser.add_argument('--overwrite', action='store_true', help='delete an existing output directory first')
    parser.add_argument('--timeout', type=float, default=120, help='per invocation seconds')
    parser.add_argument('--case', action='append', default=[], help='select case name (repeatable)')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    common_arguments(parser, 'exact-check')
    parser.add_argument('--suite', choices=('smoke', 'full', 'long', 'all'), default='full',
                        help='full includes smoke; long runs only larger fixtures; all combines both')
    parser.add_argument('--large-fixture', type=Path, help='optional supplied FASTA, exercised with -l10 (may need several GB RAM)')
    args = parser.parse_args()
    if args.timeout <= 0:
        parser.error('--timeout must be positive')
    cases = [case for case in exact_cases(args.large_fixture)
             if args.suite == 'all' or case.group == args.suite or (args.suite == 'full' and case.group == 'smoke')]
    if args.case:
        cases = [case for case in cases if case.name in args.case]
        if set(args.case) != {case.name for case in cases}:
            parser.error('unknown case or case excluded by --suite')
    reference, candidate = args.reference.resolve(), args.candidate.resolve()
    if not reference.is_file():
        print('SKIP: original TRF not found at %s; run tools/oracle/build_original_trf.sh' % reference)
        return SKIP_EXIT
    report = provenance(reference, candidate)
    report.update(suite=args.suite, cases=[], golden_policy='The small_test and diff_cn goldens are checked.')
    args.output = prepare_output(args)
    for case in cases:
        base = args.output / case.name
        reference_work, candidate_work = base / 'reference', base / 'candidate'
        ref = execute(reference, case, reference_work, args.timeout)
        cand = execute(candidate, case, candidate_work, args.timeout)
        differences = compare(ref, cand, reference_work, candidate_work)
        for label, run in [('reference', ref), ('candidate', cand)]:
            if run['exit'] != case.expected_exit:
                differences.append({'field': 'unexpected_exit', 'binary': label, 'expected': case.expected_exit, 'exit': run['exit']})
        row = {'case': case.name, 'expected_exit': case.expected_exit, 'description': case.description, 'input_sha256': sha256(reference_work / 'input.fasta'),
               'reference': ref, 'candidate': cand, 'differences': differences}
        if case.golden:
            row['upstream_golden_match'] = first_difference(case.golden, reference_work / 'stdout') is None
            if not row['upstream_golden_match']:
                differences.append({'field': 'upstream_golden', 'file': str(case.golden)})
        if case.expected_period:
            periods = {int(fields[2]) for line in (reference_work / 'stdout').read_bytes().splitlines()
                       if len(fields := line.split()) > 2 and fields[0].isdigit() and fields[2].isdigit()}
            row['expected_period_observed'] = case.expected_period in periods
            if not row['expected_period_observed']:
                differences.append({'field': 'coverage', 'reason': 'planted period was not reported', 'period': case.expected_period})
        if case.require_pagination:
            row['pagination_observed'] = any(name.endswith('.2.html') for name in ref['files'])
            row['pagination_repeat_count'] = sum(
                line.split()[0].isdigit() for path in reference_work.glob('*.dat')
                for line in path.read_bytes().splitlines() if line.split())
            if not row['pagination_observed'] or row['pagination_repeat_count'] <= 120:
                differences.append({'field': 'coverage', 'reason': 'pagination fixture did not create a second HTML page'})
        row['exact_match'] = not differences
        report['cases'].append(row)
        report['passed'] = all(item['exact_match'] for item in report['cases'])
        (args.output / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
        print(case.name, 'MATCH' if row['exact_match'] else 'FAIL',
              f"reference={ref['wall_seconds']:.3f}s candidate={cand['wall_seconds']:.3f}s", flush=True)
        if differences:
            print(json.dumps(differences, indent=2), flush=True)
    print(f"{sum(row['exact_match'] for row in report['cases'])}/{len(cases)} exact matches; {args.output / 'report.json'}")
    return 0 if report.get('passed') and cases else 1


if __name__ == '__main__':
    raise SystemExit(main())
