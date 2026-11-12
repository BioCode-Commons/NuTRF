# Benchmarks

[← NuTRF README](README.md#benchmarks)

We compare **Tandem Repeats Finder 4.10.0** (upstream commit `355c1f92b24f6b2e71650eda3549d46f0316e571`) with **NuTRF**. All numbers on this page were measured.

## Data and setup

| Test setup | Details |
| :--- | :--- |
| Datasets | Whole chromosomes chr21 (46.7 Mb), chr20 (64.4 Mb), and chr1 (249.0 Mb), extracted from the reference with `samtools faidx`; no subsetting |
| Reference | `GCA_000001405.15_GRCh38_no_alt_analysis_set.fna` |
| Hardware | AMD EPYC 9634 (84 cores, AVX2), 188 GB RAM, Linux, NVMe storage |
| Build | GCC 15.2.0, `-O2 -g -fno-align-functions -fno-align-loops`, identical for both programs |
| Run settings | One process pinned to one core, `nice 5`; several runs at a time on separate cores of a shared machine; page cache not cleared |
| Replication | Three runs per configuration; the table gives the median wall time |
| Measurement | GNU `time -v` for wall time, CPU time, and peak RSS |

Neither program is multithreaded, so there is no thread-scaling measurement. CPU time equals wall time to within 0.1 s in every run.

## Results

| Input | Arguments | Original TRF | NuTRF | Speedup | Peak RSS, original → NuTRF | Output identical |
| :--- | :--- | ---: | ---: | ---: | ---: | :--- |
| chr21 | `2 5 7 80 10 50 2000 -h -ngs` | 194.2 s | 34.5 s | **5.63×** | 215 → 218 MiB | Yes |
| chr21 | `2 7 7 80 10 50 500 -h -ngs` | 109.4 s | 24.2 s | **4.52×** | 173 → 173 MiB | Yes |
| chr21 | `2 5 7 80 10 50 2000 -d -h` | 193.8 s | 34.4 s | **5.63×** | 215 → 217 MiB | Yes |
| chr21 | `2 5 7 80 10 50 2000 -l 6 -h -ngs` | 194.8 s | 34.4 s | **5.66×** | 557 → 557 MiB | Yes |
| chr20 | `2 5 7 80 10 50 2000 -h -ngs` | 858.7 s | 138.1 s | **6.22×** | 1299 → 1296 MiB | Yes |
| chr20 | `2 7 7 80 10 50 500 -h -ngs` | 451.1 s | 79.9 s | **5.65×** | 1285 → 1281 MiB | Yes |
| chr20 | `2 5 7 80 10 50 2000 -d -h` | 858.8 s | 137.8 s | **6.23×** | 1299 → 1296 MiB | Yes |
| chr20 | `2 5 7 80 10 50 2000 -l 6 -h -ngs` | 315.8 s | 57.9 s | **5.45×** | 4068 → 4064 MiB | Yes |
| chr1 | `2 5 7 80 10 50 2000 -l 6 -h -ngs` | 630.8 s | 144.9 s | **4.35×** | 5296 → 5290 MiB | Yes |

The spread between the fastest and slowest of the three runs was below 3% in every configuration. Per-run records, including output sizes and SHA-256 digests, are in [docs/benchmarks/whole_chromosome_runs.tsv](docs/benchmarks/whole_chromosome_runs.tsv); digests of the two executables are in [docs/benchmarks/binaries.sha256](docs/benchmarks/binaries.sha256).

### Short inputs

Subsets of 150/151-base whole-genome reads and a set of 3 Mb GRCh38 windows, with `2 5 7 80 10 50 2000`, five runs each.

| Input | Original TRF | NuTRF | Speedup | Output identical |
| :--- | ---: | ---: | ---: | :--- |
| Reads, SRR7205175 | 1.59 s | 0.62 s | 2.56× | Yes |
| Reads, SRR7205167 | 3.00 s | 1.00 s | 3.00× | Yes |
| Reads, ERR1955393 | 1.56 s | 0.59 s | 2.64× | Yes |
| 3 Mb GRCh38 windows | 10.87 s | 2.42 s | 4.49× | Yes |

## Output equivalence

For every run in the tables:

- standard output (with `-ngs`) or the `.dat` file (with `-d`) matched the original byte for byte;
- standard error and exit status matched;
- the list of files written matched; and
- records were in the same order. No sorting or normalization was applied before comparison.

Further checks:

- The exact-output suite ([tests/integration/exact_check.py](tests/integration/exact_check.py), `--suite all`) matched the original executable on all 61 cases, and on a 10 Mb single-record fixture run with `-l10`.
- Builds with the AVX2 narrowband rows disabled and with all vector code disabled each reproduced the original chr21 output.
- An AddressSanitizer and UndefinedBehaviorSanitizer build reported nothing on a 5 Mb chr21 window, with and without `-l 6`.

## Limits of these measurements

- **chr1 with the default `-l` has no result.** The original was stopped after two hours without finishing. The chr1 centromere array is longer than the default 2 Mb limit, and NuTRF does not change the cost of that case; it was not timed to completion either. The chr1 row above uses `-l 6`, which the original TRF recommends for human assemblies.
- Byte identity is established for this compiler, these flags, and these inputs. It is not a proof of equivalence for all inputs. In particular, the period regression keeps a signed integer overflow from the original, whose result depends on the compiler.
- The speedups need AVX2 at run time. In single runs on chr21 with the standard arguments, a build without the AVX2 narrowband rows took about 43 s and a build without any vector code about 102 s, against 194 s for the original.
- Runs shared the machine with other jobs on other cores.
- Only human sequence was tested. Other genomes, a complete assembly in one run, other compilers, and other operating systems are not covered.

## Reproducing

Prepare the original TRF from the repository history, build NuTRF, and compare one input:

```sh
tools/oracle/build_original_trf.sh   # writes build/oracle/trf and build/oracle/source.tar
cmake -S . -B build/release -DCMAKE_BUILD_TYPE=Release
cmake --build build/release -j 8

samtools faidx GRCh38.fa chr21 > chr21.fa
/usr/bin/time -v build/oracle/trf chr21.fa 2 5 7 80 10 50 2000 -h -ngs > original.out
/usr/bin/time -v build/release/trf chr21.fa 2 5 7 80 10 50 2000 -h -ngs > nutrf.out
cmp original.out nutrf.out
```

To run the exact-output suite:

```sh
python3 tests/integration/exact_check.py --suite all
```

`ctest --test-dir build/release` runs the same suite as the test `exact_output`, together with the kernel tests.
