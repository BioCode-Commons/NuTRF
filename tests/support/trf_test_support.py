"""Paths, the original-TRF oracle, and compiler helpers shared by the tests."""
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys
import tarfile
import tempfile

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / 'src'
FIXTURES = ROOT / 'tests' / 'fixtures'
DEFAULT_CANDIDATE = ROOT / 'build' / 'release' / 'trf'
SKIP_EXIT = 77  # CTest reports this exit status as "skipped"
CFLAGS = ['-O2', '-g', '-fno-align-functions', '-fno-align-loops']


def oracle_dir():
    """Directory written by tools/oracle/build_original_trf.sh."""
    return Path(os.environ.get('TRF_ORACLE_DIR', ROOT / 'build' / 'oracle')).resolve()


def require_oracle(*names):
    """Return the oracle directory, or skip the test if a needed file is absent."""
    directory = oracle_dir()
    missing = [name for name in (names or ('source.tar', 'trf'))
               if not (directory / name).is_file()]
    if missing:
        print('SKIP: original TRF not found (%s missing in %s); '
              'run tools/oracle/build_original_trf.sh' % (', '.join(missing), directory))
        sys.exit(SKIP_EXIT)
    return directory


def original_source(name):
    """Text of one file of the original TRF sources, for example 'src/tr30dat.c'."""
    with tarfile.open(oracle_dir() / 'source.tar') as archive:
        return archive.extractfile(name).read().decode()


def extract_original_sources(directory):
    """Unpack the original src/ files into directory."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    with tarfile.open(oracle_dir() / 'source.tar') as archive:
        for member in archive.getmembers():
            path = Path(member.name)
            if member.isfile() and path.parent == Path('src'):
                (directory / path.name).write_bytes(archive.extractfile(member).read())
    return directory


def cc():
    return shlex.split(os.environ.get('CC', 'cc'))


def header_guards(keep=()):
    """Include guards of the NuTRF headers, except those listed in keep.

    A harness that supplies its own copy of the program state defines these
    guards before it includes a module source file, so that the module's
    #include lines for those headers have no effect."""
    guards = []
    for path in sorted(SRC.rglob('*.h')):
        if path.relative_to(SRC).as_posix() in keep:
            continue
        match = re.search(r'^#ifndef (TRF_\w+_H)\n#define \1$', path.read_text(), re.M)
        assert match, path
        guards.append(match.group(1))
    return guards


def define_guards(keep=()):
    return ''.join('#define %s\n' % guard for guard in header_guards(keep))


def require_sanitizers(*flags):
    """Skip the test if the compiler cannot build and run a program with flags."""
    with tempfile.TemporaryDirectory(prefix='trf-sanitizer-probe-') as directory:
        source = Path(directory) / 'probe.c'
        source.write_text('int main(void) { return 0; }\n')
        binary = Path(directory) / 'probe'
        built = subprocess.run([*cc(), *flags, str(source), '-o', str(binary)],
                               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if built.returncode or subprocess.run([str(binary)]).returncode:
            print('SKIP: the C compiler cannot build with ' + ' '.join(flags))
            sys.exit(SKIP_EXIT)
