# Included from the project root; paths are root-relative.
if(BUILD_TESTING)
    find_package(Python3 COMPONENTS Interpreter REQUIRED)

    # Differential tests compare against the original TRF 4.10.0. Run
    # tools/oracle/build_original_trf.sh once to prepare it; until then those
    # tests are reported as skipped.
    set(TRF_ORACLE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/build/oracle" CACHE PATH
        "Directory with source.tar and trf of the original TRF (tools/oracle/build_original_trf.sh)")
    option(TRF_TEST_LARGE_FIXTURE
        "Also compare the 10 Mb test_seqs fixture with -l10 (about 6 GB of memory per run)" OFF)

    set(TRF_TEST_DIR ${CMAKE_CURRENT_SOURCE_DIR}/tests)
    set(TRF_TEST_WORK ${CMAKE_CURRENT_BINARY_DIR}/test-work)
    file(MAKE_DIRECTORY ${TRF_TEST_WORK})

    # trf_python_test(<name> <script> [arguments...])
    function(trf_python_test name script)
        add_test(NAME ${name}
            COMMAND ${Python3_EXECUTABLE} ${TRF_TEST_DIR}/${script} ${ARGN})
        set_tests_properties(${name} PROPERTIES
            ENVIRONMENT "TRF_ORACLE_DIR=${TRF_ORACLE_DIR};CC=${CMAKE_C_COMPILER}"
            SKIP_RETURN_CODE 77
            TIMEOUT 1200)
    endfunction()

    # Tests that need only this build.
    add_test(NAME cli_version COMMAND trf -v)
    trf_python_test(golden_outputs integration/golden_check.py --trf $<TARGET_FILE:trf>)
    trf_python_test(harness_selftest integration/harness_selftest.py
        --reference $<TARGET_FILE:trf>)
    trf_python_test(similarity unit/test_similarity.py)
    trf_python_test(narrowband_wrap_proofs unit/test_narrowband.py --wrap-only)
    trf_python_test(narrowband_simd_guards unit/test_narrowband_simd_guards.py)
    trf_python_test(narrowband_simd_rows unit/test_narrowband_simd_rows.py)
    trf_python_test(newwrap_kernel_self unit/test_newwrap.py --kernel-only)

    # Differential tests against the original TRF.
    trf_python_test(exact_output integration/exact_check.py --suite all
        --candidate $<TARGET_FILE:trf>
        --output ${TRF_TEST_WORK}/exact-check --overwrite)
    trf_python_test(periods unit/test_periods.py)
    trf_python_test(list_sort unit/test_list_sort.py)
    trf_python_test(narrowband unit/test_narrowband.py)
    trf_python_test(narrowband_state unit/test_narrowband_state.py)
    trf_python_test(narrowband_simd_avx2 unit/test_narrowband_simd.py --variant avx2)
    trf_python_test(narrowband_simd_sse_ramp unit/test_narrowband_simd.py --variant sse-ramp)
    trf_python_test(narrowband_simd_packed unit/test_narrowband_simd.py --variant packed)
    trf_python_test(narrowband_simd_force_scalar unit/test_narrowband_simd.py --force-scalar)
    trf_python_test(narrowband_simd_portable unit/test_narrowband_simd.py --portable)
    trf_python_test(newwrap unit/test_newwrap.py)
    set_tests_properties(exact_output periods list_sort narrowband narrowband_state
        narrowband_simd_avx2 narrowband_simd_sse_ramp narrowband_simd_packed
        narrowband_simd_force_scalar narrowband_simd_portable newwrap
        PROPERTIES LABELS oracle)

    if(TRF_TEST_LARGE_FIXTURE)
        add_test(NAME large_fixture_generate
            COMMAND ${Python3_EXECUTABLE} ${TRF_TEST_DIR}/support/make_test_seqs.py
                ${TRF_TEST_WORK}/test_seqs.fasta)
        set_tests_properties(large_fixture_generate PROPERTIES FIXTURES_SETUP trf_large_fixture)
        trf_python_test(exact_output_large integration/exact_check.py --suite long
            --candidate $<TARGET_FILE:trf>
            --large-fixture ${TRF_TEST_WORK}/test_seqs.fasta --timeout 1800
            --output ${TRF_TEST_WORK}/exact-check-large --overwrite)
        set_tests_properties(exact_output_large PROPERTIES
            FIXTURES_REQUIRED trf_large_fixture LABELS "oracle;large" TIMEOUT 3600
            RUN_SERIAL TRUE)
    endif()
endif()
