#!/usr/bin/env python3
"""Prove that real candidate processes with altered outputs/status cannot pass."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'support'))
from trf_test_support import DEFAULT_CANDIDATE


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference', type=Path, default=DEFAULT_CANDIDATE,
                        help='any working TRF executable; the faults are injected around it')
    args = parser.parse_args()
    reference = args.reference.resolve()
    injections = {
        'stdout': ("sys.stdout.buffer.write(b'INJECTED\\n')", 'file_bytes'),
        'stderr': ("sys.stderr.buffer.write(b'INJECTED\\n')", 'file_bytes'),
        'filename': ("Path('injected.dat').write_bytes(b'INJECTED')", 'filenames'),
        'generated-bytes': ("next(Path('.').glob('*.dat')).open('ab').write(b'INJECTED')", 'file_bytes'),
        'exit': ('status = 23', 'exit'),
        'timeout': ('time.sleep(10)', 'incomplete'),
        'identical-signal': ('os.kill(os.getpid(), signal.SIGTERM)', 'signal_exit'),
    }
    with tempfile.TemporaryDirectory(prefix='trf-harness-selftest-') as tmp:
        root = Path(tmp)
        for name, (mutation, expected) in injections.items():
            candidate = root / (name + '.py')
            # Explicit executable= preserves argv[0] inside the real TRF child.
            candidate.write_text(f'#!{sys.executable}\n'
                                 'from pathlib import Path\nimport subprocess, sys, time, os, signal\n'
                                 f"status = subprocess.run(['trf', *sys.argv[1:]], executable={str(reference)!r}).returncode\n"
                                 + mutation + '\nraise SystemExit(status)\n')
            candidate.chmod(0o755)
            output = root / name
            command = [sys.executable, str(Path(__file__).resolve().with_name('exact_check.py')), '--reference', str(candidate if name == 'identical-signal' else reference),
                       '--candidate', str(candidate), '--case', 'small_test-data', '--output', str(output),
                       '--timeout', '0.5' if name == 'timeout' else '5']
            result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=15)
            report = json.loads((output / 'report.json').read_text())
            differences = report['cases'][0]['differences']
            if result.returncode != 1 or expected not in {item['field'] for item in differences}:
                print(name, 'FAILED', result.stdout.decode(), result.stderr.decode(), differences)
                return 1
            print(name, 'correctly rejected:', ', '.join(sorted({item['field'] for item in differences})), flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
