"""Deterministic synthetic inputs shared by exact checks and benchmarks."""
from dataclasses import dataclass
from pathlib import Path
import random
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'support'))
from trf_test_support import FIXTURES, ROOT

SEED = 20261104
DEFAULT_PARAMS = ('2', '5', '7', '80', '10', '50', '2000')


@dataclass
class Case:
    name: str
    data: bytes = b''
    flags: tuple = ('-ngs',)
    params: tuple = DEFAULT_PARAMS
    stdin: bool = False
    argv: tuple | None = None
    group: str = 'full'
    description: str = ''
    golden: Path | None = None
    require_pagination: bool = False
    expected_period: int | None = None
    expected_exit: int = 0
    source: Path | None = None

    def arguments(self):
        return list(self.argv) if self.argv is not None else [
            '-' if self.stdin else 'input.fasta', *self.params, *self.flags]


def dna(rng, length):
    return ''.join(rng.choices('ACGT', k=length))


def fasta(records):
    return ''.join('>' + name + '\n' + seq + '\n' for name, seq in records).encode('ascii')


def mutate(rng, seq, substitutions=29, indels=0):
    result = []
    for i, base in enumerate(seq):
        if indels and i % indels == indels - 1:
            if rng.randrange(2):
                result.append(rng.choice('ACGT'))
            else:
                continue
        result.append(rng.choice('ACGT'.replace(base, '')) if i % substitutions == 0 else base)
    return ''.join(result)


def pagination_input():
    rng = random.Random(SEED + 120)
    # Independent motifs separated by ambiguous bases make >120 distinct hits.
    return fasta([('pagination', ('N' * 100).join(dna(rng, 17) * 9 for _ in range(135)))])


def exact_cases(large_fixture=None):
    rng = random.Random(SEED)
    cases = []
    for name in ('small_test', 'diff_cn'):
        data = (FIXTURES / (name + '.fasta')).read_bytes()
        for mode, flags in [('ngs', ('-ngs',)), ('data', ('-h',)),
                            ('all', ('-d', '-m', '-f')), ('raw', ('-ngs', '-r'))]:
            cases.append(Case(name + '-' + mode, data, flags, group='smoke',
                              golden=FIXTURES / (name + '.fasta.out') if mode == 'ngs' else None))
    small = (FIXTURES / 'small_test.fasta').read_bytes()
    cases += [Case('stdin', small, stdin=True, group='smoke'),
              Case('no-repeat', fasta([('none', dna(rng, 8000))]), group='smoke'),
              Case('lowercase-ambiguous', fasta([('lower', ('acgttgca' * 30) + 'nryswkmbdhvNX!?0123' + ('ta' * 80))])),
              Case('multiple-records', fasta([('empty', ''), ('repeat', 'ACTG' * 60),
                                            ('random', dna(rng, 700)), ('ambiguous', 'N' * 150)])),
              Case('empty-input', b'', group='smoke', expected_exit=1),
              Case('empty-record', b'>empty\n'),
              Case('malformed-no-header', b'ACGTACGTACGT\n', expected_exit=1),
              Case('malformed-headers', b'>one\r\nACGTACGT\r\n>two\n>three\nNNN\n'),
              Case('version', argv=('-v',), group='smoke'),
              Case('usage', argv=(), group='smoke', expected_exit=1),
              Case('stdin-without-ngs', small, flags=('-h',), stdin=True, expected_exit=255),
              Case('missing-input', argv=('missing.fasta', *DEFAULT_PARAMS, '-ngs'), expected_exit=1)]
    for label, index, value in [('period-zero', 6, '0'), ('period-over', 6, '2001'),
                                ('match-invalid', 0, 'abc'), ('match-negative', 0, '-1'),
                                ('score-invalid', 5, '999999999999999999999')]:
        params = list(DEFAULT_PARAMS)
        params[index] = value
        cases.append(Case(label, small, params=tuple(params), expected_exit=1))
    cases.append(Case('invalid-option', small, flags=('-ngs', '--unknown-option')))
    for period in (1, 2, 20, 21, 50, 51, 100, 101, 499, 500, 501, 1999, 2000):
        motif = dna(rng, period)
        data = fasta([('period-' + str(period), dna(rng, 70) + motif * max(5, 120 // period) + dna(rng, 70))])
        params = (*DEFAULT_PARAMS[:6], str(period))
        cases.append(Case('period-' + str(period), data, params=params,
                          description='Planted repeat and MaxPeriod at boundary ' + str(period), expected_period=period))
    # Keep the same mixed input while varying MaxPeriod to exercise cutoff decisions.
    mixed = fasta([('mixed-boundaries', ('N' * 100).join(
        dna(rng, p) * 5 for p in (20, 21, 50, 51, 100, 101, 499, 500, 501, 1999, 2000)))])
    for maximum in (20, 21, 50, 51, 100, 101, 500, 2000):
        cases.append(Case('maxperiod-' + str(maximum), mixed,
                          params=(*DEFAULT_PARAMS[:6], str(maximum))))
    cases += [Case('pm75-pi20-mixed-baseline-failure', mixed, params=('2', '5', '7', '75', '20', '50', '2000'), expected_exit=255, description='Input-specific inherited exit255 with empty output on long mixed periods and PI20'),
              Case('pm75-pi20', small, params=('2', '5', '7', '75', '20', '50', '2000')),
              Case('alternate-scoring', small, params=('2', '7', '7', '80', '10', '50', '2000')),
              Case('crlf', small.replace(b'\n', b'\r\n')),
              Case('no-final-newline', small.rstrip(b'\r\n')),
              Case('overlap', fasta([('overlap', 'ACGT' * 80 + 'ACGTACGA' * 60 + 'ACGT' * 80)])),
              Case('pagination-html', pagination_input(), flags=('-d', '-m', '-f'), require_pagination=True)]
    for seed in range(6):
        local = random.Random(SEED + 1000 + seed)
        period = (7, 20, 51, 101, 500, 2000)[seed]
        seq = mutate(local, dna(local, period) * max(6, 6000 // period),
                     substitutions=17 + seed, indels=31 + seed)
        cases.append(Case('substitution-indel-' + str(seed),
                          fasta([('seed-' + str(SEED + 1000 + seed), seq)])))
    random_seq = dna(rng, 200000)
    imperfect = mutate(rng, dna(rng, 125) * 1600, substitutions=23)
    cases.append(Case('synthetic-long', fasta([('random', random_seq), ('imperfect', imperfect)]), group='long'))
    # The source parses av[8] instead of optarg for -l, as the original does.
    if large_fixture:
        cases.append(Case('supplied-10mb-l10', flags=('-ngs', '-l10'),
                          source=Path(large_fixture).resolve(),
                          group='long', description='Supplied fixture with -l10 and MaxPeriod=2000; -l parser quirk of the original retained'))
    return cases


def original_mixed():
    rng = random.Random(SEED)
    random_seq = dna(rng, 200000)
    motif = dna(rng, 125)
    repeated = list(motif * 1600)
    for i in range(0, len(repeated), 23):
        repeated[i] = rng.choice('ACGT')
    return fasta([('random', random_seq), ('imperfect', ''.join(repeated))])


def benchmark_cases():
    rng = random.Random(SEED + 7000)
    return [
        Case('original-mixed-400k', original_mixed()),
        Case('indel-rich-ngs', fasta([('synthetic-indel-rich', mutate(rng, dna(rng, 125) * 1200, substitutions=29, indels=31))])),
        Case('random-ngs', fasta([('synthetic-random', dna(rng, 400000))])),
        Case('repeat-ngs', fasta([('synthetic-imperfect', mutate(rng, dna(rng, 125) * 1600, substitutions=23))])),
        Case('short-records-ngs', fasta([(f'synthetic-{i}', dna(rng, 31) * 8) for i in range(300)])),
        Case('output-data', pagination_input(), flags=('-h', '-m', '-f')),
        Case('output-html', pagination_input(), flags=('-d', '-m', '-f')),
    ]
