#!/usr/bin/env python3
"""Run TRF on the bundled FASTA fixtures and compare standard output with the
expected -ngs output stored next to each fixture."""
import argparse
from pathlib import Path
import subprocess
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'support'))
from trf_test_support import DEFAULT_CANDIDATE, FIXTURES

ARGUMENTS = ('2', '5', '7', '80', '10', '50', '2000', '-ngs')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--trf', type=Path, default=DEFAULT_CANDIDATE)
    args = parser.parse_args()
    failed = False
    for name in ('small_test', 'diff_cn'):
        fixture = FIXTURES / (name + '.fasta')
        expected = (FIXTURES / (name + '.fasta.out')).read_bytes()
        with tempfile.TemporaryDirectory(prefix='trf-golden-') as work:
            result = subprocess.run([str(args.trf.resolve()), str(fixture), *ARGUMENTS],
                                    cwd=work, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        matched = result.returncode == 0 and result.stdout == expected
        print(name, 'MATCH' if matched else 'FAIL (exit %d)' % result.returncode)
        failed = failed or not matched
    return 1 if failed else 0


if __name__ == '__main__':
    raise SystemExit(main())
