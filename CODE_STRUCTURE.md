# Reading the NuTRF code

[← NuTRF README](README.md#implementation-status-and-roadmap)

Start with [src/cli/main.c](src/cli/main.c). It parses the command line into a parameter set and calls `TRFControlRoutine`. Everything else is reached from there.

The source is one directory per responsibility under [src/](src/). Each `.c` file has a header of the same name that declares what other modules may call. Includes are written relative to `src/`, for example `#include "alignment/newwrap.h"`.

| Directory | Responsibility |
| --- | --- |
| [src/cli/](src/cli/) | `main`, argument parsing, usage text and banner |
| [src/run/](src/run/) | Control of a run and of one sequence; FASTA loading |
| [src/core/](src/core/) | Shared types, constants, and program state; messages; work-array allocation |
| [src/detection/](src/detection/) | k-tuple scan, distance lists, statistical criteria, criteria tables |
| [src/periods/](src/periods/) | Period ranking and the list of best periods |
| [src/alignment/](src/alignment/) | Similarity matrix, small-period and narrowband wraparound alignment, traceback |
| [src/consensus/](src/consensus/) | Consensus pattern and repeat statistics |
| [src/redundancy/](src/redundancy/) | List sorts and redundancy elimination |
| [src/output/](src/output/) | Alignment text, HTML, `.dat`, `-ngs`, and mask writers |

## From a sequence to a repeat table

```text
command line -> load FASTA record -> k-tuple scan -> statistical criteria
    -> rank candidate periods -> wraparound alignment (small or narrowband)
    -> consensus and statistics -> sort and remove redundancy -> write output
```

| Step | Read these files | What the code does |
| --- | --- | --- |
| Parse the command | [src/cli/main.c](src/cli/main.c), [src/cli/arguments.c](src/cli/arguments.c), [src/cli/usage.c](src/cli/usage.c) | `ParseArguments` reads the seven numeric arguments and the options, and fills the global parameter set. `main` then runs the control routine and reports the status. |
| Drive a run | [src/run/control.c](src/run/control.c), [src/run/fasta.c](src/run/fasta.c) | `TRFControlRoutine` loops over FASTA records and opens the data, mask, and summary files. `LoadSequenceFromFileEugene` reads one record. |
| Process a sequence | [src/run/sequence.c](src/run/sequence.c) | `TRF` allocates the work arrays, runs the scan, and resolves the output of one sequence. |
| Scan for candidates | [src/detection/tuple_scan.c](src/detection/tuple_scan.c), [src/detection/distance_list.c](src/detection/distance_list.c), [src/detection/distance_seen.c](src/detection/distance_seen.c) | `newtupbo` slides k-tuples of several sizes along the sequence and records distances between equal tuples. It calls the alignment, consensus, and statistics steps for each accepted candidate. |
| Apply the criteria | [src/detection/criteria.c](src/detection/criteria.c), [src/detection/coin_toss.c](src/detection/coin_toss.c) | `new_meet_criteria_3` and `multiples_criteria_4` decide whether a distance is aligned and reported. `coin_toss.c` holds the waiting-time and k-run-sum tables. |
| Rank periods | [src/periods/period_ranking.c](src/periods/period_ranking.c), [src/periods/best_period_list.c](src/periods/best_period_list.c) | `GetTopPeriods` in `criteria.c` forwards to `trf_get_top_periods`, a function without global state that counts distances between equal dinucleotides and ranks the most likely period sizes. |
| Align, small periods | [src/alignment/newwrap.c](src/alignment/newwrap.c) | `newwrap` aligns the sequence against a pattern with wraparound dynamic programming over the whole pattern width. |
| Align, large periods | [src/alignment/narrowband.c](src/alignment/narrowband.c), [src/alignment/narrowband_scalar.c](src/alignment/narrowband_scalar.c), [src/alignment/narrowband_simd.c](src/alignment/narrowband_simd.c) | `narrowbandwrap` fills only a band around the expected diagonal. `narrowband.c` selects the kernel for a call; the other two files hold the scalar kernel and the SSE4.1 and AVX2 row kernels. |
| Trace back | [src/alignment/pair_alignment.c](src/alignment/pair_alignment.c), [src/alignment/similarity.c](src/alignment/similarity.c) | `get_pair_alignment_with_copynumber` and its narrowband variant trace the alignment back into `AlignPair` and compute the copy number. `init_sm` builds the similarity matrix. |
| Describe a repeat | [src/consensus/consensus.c](src/consensus/consensus.c), [src/consensus/statistics.c](src/consensus/statistics.c) | `get_consensus` and `get_statistics` compute the consensus and the percentages, and append the repeat to the repeat list. |
| Remove redundancy | [src/redundancy/clean.c](src/redundancy/clean.c), [src/redundancy/list_sort.c](src/redundancy/list_sort.c), [src/redundancy/redundancy.c](src/redundancy/redundancy.c) | `TRFClean` sorts the repeat list and drops repeats made redundant by others. The list sorts are stable merge sorts. |
| Write results | [src/output/alignment_text.c](src/output/alignment_text.c), [src/output/html.c](src/output/html.c), [src/output/summary_html.c](src/output/summary_html.c), [src/output/data_file.c](src/output/data_file.c), [src/output/mask.c](src/output/mask.c) | Writes the alignment text, HTML tables, the summary page, the `.dat` file or `-ngs` stream, and the masked sequence. |

Shared types and constants are in [src/core/types.h](src/core/types.h) and [src/core/constants.h](src/core/constants.h). The program state is declared in [src/core/globals.h](src/core/globals.h) and defined once in [src/core/globals.c](src/core/globals.c). The similarity matrix `SM` is written only by `init_sm`; the alignment kernels rely on that.

## What NuTRF changed, and where

The original TRF is a single translation unit: `trf.c` includes `trfrun.h`, `tr30dat.h`, `tr30dat.c`, and `trfclean.h`. NuTRF compiles the same routines as the modules above. The table lists the files whose code differs from the original; the other files hold the original's code.

| File | Original counterpart | Difference |
| --- | --- | --- |
| [src/periods/period_ranking.c](src/periods/period_ranking.c) | `GetTopPeriods` in `tr30dat.c` | Bucketed positions with a moving lower bound; AVX2 byte-compare histogram for patterns of at least 128 bases. |
| [src/alignment/narrowband_scalar.c](src/alignment/narrowband_scalar.c), [narrowband_scalar.h](src/alignment/narrowband_scalar.h) | `narrowbandwrap` and its macros in `tr30dat.c` | Conditional wrap instead of `%`, and direct reads of the previous row. |
| [src/alignment/narrowband_simd.c](src/alignment/narrowband_simd.c), [src/alignment/narrowband.c](src/alignment/narrowband.c) | None | Vector row kernels: prefix-maximum form of the left-gap recurrence, AVX2 loops for runs of ordinary rows, per-row score-range check; selection of the kernel at run time. |
| [src/alignment/newwrap.c](src/alignment/newwrap.c) | `newwrap` in `tr30dat.c` | Reuses the first circular pass, skips the preliminary forward scan, and has AVX2 rows. The original function is kept in the file as the path for unsupported scoring. |
| [src/alignment/similarity.c](src/alignment/similarity.c) | `init_sm` in `tr30dat.c` | Rebuilds the matrix only when the scoring parameters change. |
| [src/detection/tuple_scan.c](src/detection/tuple_scan.c), [src/detection/distance_list.c](src/detection/distance_list.c) | `newtupbo` and the distance-list routines in `tr30dat.c` | Division-free ring indices and tuple codes; table lookup for the ACGT test. |
| [src/redundancy/list_sort.c](src/redundancy/list_sort.c) | `SortByIndex` and `SortByCount` in `trfclean.h` | Stable merge sort in place of bubble sorts. |
| [src/output/data_file.c](src/output/data_file.c), [src/output/summary_html.c](src/output/summary_html.c) | Inline code of `TRFControlRoutine` in `trfrun.h` | The same statements as functions, called for single-sequence and multi-sequence input. |

Each vector path has a scalar path that produces the same integers. The CPU is tested at run time; these compile-time switches force the slower paths for testing, for example with `-DCMAKE_C_FLAGS=-DTRF_NARROWBAND_DISABLE_AVX2`:

| Switch | Effect |
| --- | --- |
| `-DTRF_NARROWBAND_DISABLE_AVX2` | Use the SSE4.1 narrowband rows |
| `-DTRF_NARROWBAND_DISABLE_RAMP` | With the switch above, use the packed SSE4.1 scan |
| `-DTRF_NARROWBAND_DISABLE_SIMD`, `-DTRF_NARROWBAND_FORCE_SCALAR` | Use the scalar narrowband kernel |
| `-DTRF_NEWWRAP_DISABLE_SIMD`, `-DTRF_NEWWRAP_FORCE_SCALAR` | Use the scalar small-period rows |
| `-DTRF_PERIODS_DISABLE_SIMD`, `-DTRF_PERIODS_FORCE_SCALAR` | Use the scalar period count |

## What each kernel preserves

Output identity depends on properties of the original code that are easy to break. The kernels keep all of them.

| Kernel | Preserved from the original |
| --- | --- |
| Period ranking | Every raw distance count; the rule that the nearest predecessor at or beyond 6,000 bases still counts, so a sparse gap can contribute a larger distance; the regression expressions and their accumulation order; strict `>` scans, so equal peaks keep the smallest period first; the signed overflow in the regression. |
| Narrowband alignment | Cell values, cell order within the band, recentering of the band, sentinels, endpoints, and traceback. Ties between maxima resolve to the same cell. Cells outside the band are left as the original leaves them. |
| Small-period alignment | The alignment endpoint, the live matrix region, and the traceback. Scratch arrays hold the same values on exit. |
| Redundancy elimination | The order of equal elements, because the sorts are stable. |

## Build

[CMakeLists.txt](CMakeLists.txt) builds one static library per source directory and links them into `trf`. The Release flags are `-O2 -g -fno-align-functions -fno-align-loops`, the flags of the original TRF build.

The k-tuple scan calls small functions of other source files in the same directory for every tuple match. With `TRF_UNITY_BUILD=ON`, the default, CMake compiles the sources of each directory as one translation unit, so the compiler can inline those calls. `TRF_ENABLE_RELEASE_IPO=ON` enables link-time optimization across all modules instead of, or in addition to, the unity build. The exact-output suite matches the original in each of these configurations.

## Tests

Tests are Python scripts under [tests/](tests/) that use only the standard library and a C compiler. [tests/CMakeLists.txt](tests/CMakeLists.txt) registers them with CTest. Most compare against the original TRF. Run [tools/oracle/build_original_trf.sh](tools/oracle/build_original_trf.sh) once to create `build/oracle/source.tar` and `build/oracle/trf` from the upstream commit in the repository history; CTest passes that directory to the scripts (`TRF_ORACLE_DIR`). Without it, the differential tests exit with status 77 and CTest reports them as skipped.

```sh
tools/oracle/build_original_trf.sh
ctest --test-dir build/release --output-on-failure
```

| Script | CTest name | What it checks |
| --- | --- | --- |
| [tests/integration/exact_check.py](tests/integration/exact_check.py), [trf_corpus.py](tests/integration/trf_corpus.py) | `exact_output` | Runs the original executable and NuTRF on a corpus of inputs and options; compares exit status, both streams, file names, and every output byte. |
| [tests/integration/golden_check.py](tests/integration/golden_check.py) | `golden_outputs` | The bundled fixtures against their stored `-ngs` output. Needs no original TRF. |
| [tests/integration/harness_selftest.py](tests/integration/harness_selftest.py) | `harness_selftest` | Injects faults to confirm that the exact-output comparison detects differences. Needs no original TRF. |
| [tests/unit/test_periods.py](tests/unit/test_periods.py) | `periods` | Period histograms, regression slope, and rankings against the original function. |
| [tests/unit/test_narrowband.py](tests/unit/test_narrowband.py), [test_narrowband_state.py](tests/unit/test_narrowband_state.py) | `narrowband`, `narrowband_state`, `narrowband_wrap_proofs` | The scalar narrowband kernel: source text against the original function, complete matrix and global state, and end-to-end output of the original program with only this kernel replaced. The wrap proofs need no original TRF. |
| [tests/unit/test_narrowband_simd.py](tests/unit/test_narrowband_simd.py) | `narrowband_simd_avx2`, `_sse_ramp`, `_packed`, `_force_scalar`, `_portable` | The vector narrowband kernels, each variant: complete matrix and global state against the original function, including out-of-band cells. |
| [tests/unit/test_narrowband_simd_rows.py](tests/unit/test_narrowband_simd_rows.py), [test_narrowband_simd_guards.py](tests/unit/test_narrowband_simd_guards.py) | `narrowband_simd_rows`, `narrowband_simd_guards` | Single vector rows against an ordered scalar recurrence, and the boundaries of the kernel-selection test, under sanitizers. Need no original TRF. |
| [tests/unit/test_newwrap.py](tests/unit/test_newwrap.py) | `newwrap`, `newwrap_kernel_self` | Small-period alignment: every row of every scan, endpoints, live matrix, and traceback; end-to-end output of the original program with only this kernel replaced. `newwrap_kernel_self` uses the copy of the original function in `newwrap.c` and needs no original TRF. |
| [tests/unit/test_similarity.py](tests/unit/test_similarity.py) | `similarity` | The similarity matrix after repeated and changed scoring parameters. |
| [tests/unit/test_list_sort.py](tests/unit/test_list_sort.py) | `list_sort` | The merge sorts against the original bubble sorts, including stability. |

The tests reach the code in three ways. `test_similarity.py`, `test_list_sort.py`, `test_periods.py`, and `test_narrowband_simd_guards.py` compile a small harness and link it with the module's source file. The narrowband and small-period state tests compare every byte of the program state, so their harness holds that state in one structure and compiles the complete module source against it, with the module's includes of the shared headers disabled through their include guards ([tests/support/trf_test_support.py](tests/support/trf_test_support.py)); the row test reaches the static row helpers of `narrowband_simd.c` the same way. The state tests also insert observation counters into the module text to confirm that each kernel was exercised. The original function always comes from `build/oracle/source.tar`.

The exact-output corpus has 60 regular cases and a 400 kb synthetic input. It covers both small fixtures in every output mode, standard input, lowercase and ambiguous symbols, multiple and empty records, malformed input, and argument errors. A larger FASTA can be added with `--large-fixture`; `-DTRF_TEST_LARGE_FIXTURE=ON` registers that run with the 10 Mb fixture written by [tests/support/make_test_seqs.py](tests/support/make_test_seqs.py).

## Other files

| Path | Purpose |
| --- | --- |
| [CMakeLists.txt](CMakeLists.txt), [tests/CMakeLists.txt](tests/CMakeLists.txt) | Build, install, and test registration. |
| [tools/oracle/build_original_trf.sh](tools/oracle/build_original_trf.sh) | Build the original TRF for the differential tests. |
| [tools/benchmark/benchmark.py](tools/benchmark/benchmark.py) | Paired timing runs on synthetic inputs that verify output before accepting a timing. |
| [tests/fixtures/](tests/fixtures/) | Small FASTA fixtures with their expected `-ngs` output, and the table from which `make_test_seqs.py` writes `test_seqs.fasta`, the example input of the original manual. |
| [docs/TRF_README.md](docs/TRF_README.md) | The original manual: output explanations, the detection model, and the version history. |
| [docs/benchmarks/](docs/benchmarks/) | Whole-chromosome run records. |
| [docs/upstream/](docs/upstream/) | `NEWS`, `ChangeLog`, `history.txt`, and `INSTALL` of the original TRF, unchanged. |
| `AUTHORS`, `COPYING` | Upstream authorship and license. |
