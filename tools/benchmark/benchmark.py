#!/usr/bin/env python3
"""Alternating TRF trials with verified outputs and per-child wall/CPU/RSS."""
import argparse
import json
from pathlib import Path
import statistics
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'tests' / 'integration'))
from exact_check import common_arguments, compare, execute, prepare_output, provenance, sha256
from trf_corpus import benchmark_cases


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    common_arguments(parser, 'benchmark')
    parser.add_argument('--trials', type=int, default=5)
    parser.add_argument('--warmups', type=int, default=1, help='verified runs of each binary, excluded from medians')
    args = parser.parse_args()
    if args.trials < 1 or args.warmups < 0 or args.timeout <= 0:
        parser.error('trials and timeout must be positive; warmups cannot be negative')
    cases = benchmark_cases()
    if args.case:
        cases = [case for case in cases if case.name in args.case]
        if set(args.case) != {case.name for case in cases}:
            parser.error('unknown benchmark case')
    reference, candidate = args.reference.resolve(), args.candidate.resolve()
    report = provenance(reference, candidate)
    report.update(trials=args.trials, warmups=args.warmups, cases=[], passed=True,
                  corpus_kind='Deterministically generated synthetic DNA; not a real genome',
                  measurement='wall perf_counter around spawn/wait; wait4 child user/system CPU and Linux ru_maxrss KiB; output retained and verified after every pair',
                  ordering='Reference first on even numbered pairs, candidate first on odd numbered pairs; warmups included in ordering')
    args.output = prepare_output(args)
    for case in cases:
        row = {'case': case.name, 'runs': [], 'input_sha256': None}
        report['cases'].append(row)
        canonical, canonical_work = None, None
        for trial in range(args.warmups + args.trials):
            phase = 'warmup' if trial < args.warmups else 'measured'
            order = [('reference', reference), ('candidate', candidate)]
            if trial % 2:
                order.reverse()
            runs, works = {}, {}
            for label, binary in order:
                work = args.output / case.name / f'{trial:02d}-{phase}' / label
                runs[label] = execute(binary, case, work, args.timeout)
                works[label] = work
            if canonical is None:
                canonical, canonical_work = runs['reference'], works['reference']
                row['input_sha256'] = sha256(canonical_work / 'input.fasta')
            differences = compare(runs['reference'], runs['candidate'], works['reference'], works['candidate'])
            # Repeated runs must also match the first reference run, not merely each other.
            for label in ('reference', 'candidate'):
                repeated = compare(canonical, runs[label], canonical_work, works[label])
                if repeated:
                    differences.append({'field': 'repeatability', 'binary': label, 'differences': repeated})
                if runs[label]['exit'] != 0:
                    differences.append({'field': 'unexpected_exit', 'binary': label, 'exit': runs[label]['exit']})
            row['runs'].append({'trial': trial, 'phase': phase, 'order': [label for label, _ in order],
                                'exact_match': not differences, 'differences': differences, **runs})
            report['passed'] = report['passed'] and not differences
            (args.output / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
            print(case.name, phase, trial, 'MATCH' if not differences else 'FAIL',
                  ' '.join(f"{label}={runs[label]['wall_seconds']:.4f}s" for label in ('reference', 'candidate')), flush=True)
            if differences:
                print(json.dumps(differences, indent=2), flush=True)
                print('Stopping benchmark: output verification failed; partial report retained.', flush=True)
                return 1
        row['medians'] = {}
        for label in ('reference', 'candidate'):
            measured = [run[label] for run in row['runs'] if run['phase'] == 'measured']
            row['medians'][label] = {
                key: statistics.median(run[key] for run in measured)
                for key in ('wall_seconds', 'cpu_user_seconds', 'cpu_system_seconds', 'peak_rss_kib')}
            row['medians'][label]['cpu_seconds'] = statistics.median(
                run['cpu_user_seconds'] + run['cpu_system_seconds'] for run in measured)
        row['wall_speedup_reference_over_candidate'] = row['medians']['reference']['wall_seconds'] / row['medians']['candidate']['wall_seconds']
        print(case.name, 'medians', json.dumps(row['medians'], sort_keys=True), flush=True)
        (args.output / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    return 0 if report['passed'] and cases else 1


if __name__ == '__main__':
    raise SystemExit(main())
