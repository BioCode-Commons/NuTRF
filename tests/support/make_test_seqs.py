#!/usr/bin/env python3
"""Write the test_seqs.fasta fixture described by tests/fixtures/test_seqs.tsv.

Each row of the table gives a sequence identifier, a pattern, and a number of
exact tandem copies. The last record is a 10 Mb array of a 125-base pattern,
which needs -l 10. The original TRF manual (docs/TRF_README.md) uses this file
as its example input."""
import argparse
import hashlib
from pathlib import Path

TABLE = Path(__file__).resolve().parents[1] / 'fixtures' / 'test_seqs.tsv'
SHA256 = '8e784706f40f9d3924a33616c3b2a649378df287e516d3caec2dc47ac1715384'
WIDTH = 70


def records():
    lines = TABLE.read_text().splitlines()
    assert lines[0].split('\t') == ['SID', 'Pattern', 'NumCopies']
    for line in lines[1:]:
        name, pattern, copies = line.split('\t')
        yield name, pattern * int(copies)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path, help='FASTA file to write')
    args = parser.parse_args()
    digest = hashlib.sha256()
    with args.output.open('wb') as handle:
        for name, sequence in records():
            text = '>' + name + '\n' + ''.join(
                sequence[i:i + WIDTH] + '\n' for i in range(0, len(sequence), WIDTH))
            handle.write(text.encode('ascii'))
            digest.update(text.encode('ascii'))
    if digest.hexdigest() != SHA256:
        raise SystemExit('unexpected content: ' + digest.hexdigest())
    print(args.output)


if __name__ == '__main__':
    main()
