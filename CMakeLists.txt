cmake_minimum_required(VERSION 3.16)
project(NuTRF VERSION 4.10.0 LANGUAGES C)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Output identity with the original TRF was verified at -O2. Do not raise the
# optimization level or add -march or fast-math flags if identical output is
# required; SSE4.1 and AVX2 code is compiled per function and selected at run
# time, so no architecture flag is needed.
set(TRF_RELEASE_FLAGS "-O2 -g" CACHE STRING "Compiler flags of the Release build")
set(CMAKE_C_FLAGS_RELEASE "${TRF_RELEASE_FLAGS}" CACHE STRING
    "Set from TRF_RELEASE_FLAGS" FORCE)

# The k-tuple scan calls small functions of other source files for every
# tuple match. A unity build compiles the sources of each module directory as
# one translation unit so that those calls are inlined.
option(TRF_UNITY_BUILD "Compile each module directory as one translation unit" ON)

# Alternative to the unity build: link-time inlining across all modules.
option(TRF_ENABLE_RELEASE_IPO "Enable interprocedural optimization in Release" OFF)
if(TRF_ENABLE_RELEASE_IPO)
    include(CheckIPOSupported)
    check_ipo_supported(RESULT TRF_IPO_SUPPORTED OUTPUT TRF_IPO_ERROR LANGUAGES C)
    if(TRF_IPO_SUPPORTED)
        set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
    else()
        message(STATUS "Release IPO unavailable: ${TRF_IPO_ERROR}")
    endif()
endif()

include(CTest)
include(CheckCCompilerFlag)
include(CheckLibraryExists)

# Shared settings: include root, version string, and the alignment flags that
# the original TRF build uses.
add_library(trf_settings INTERFACE)
target_include_directories(trf_settings INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_compile_definitions(trf_settings INTERFACE TRF_VERSION_STRING="${PROJECT_VERSION}")
foreach(flag -fno-align-functions -fno-align-loops)
    string(MAKE_C_IDENTIFIER "TRF_HAVE${flag}" flag_variable)
    check_c_compiler_flag(${flag} ${flag_variable})
    if(${flag_variable})
        target_compile_options(trf_settings INTERFACE ${flag})
        target_link_options(trf_settings INTERFACE ${flag})
    endif()
endforeach()
check_library_exists(m floor "" TRF_HAVE_LIBM)
if(TRF_HAVE_LIBM)
    target_link_libraries(trf_settings INTERFACE m)
endif()

# One library per source directory.
function(trf_module name)
    add_library(trf_${name} STATIC ${ARGN})
    target_link_libraries(trf_${name} PUBLIC trf_settings)
    set_target_properties(trf_${name} PROPERTIES
        UNITY_BUILD ${TRF_UNITY_BUILD} UNITY_BUILD_BATCH_SIZE 0)
endfunction()

trf_module(core
    src/core/globals.c src/core/report.c src/core/arrays.c src/core/repeat_list.c)

trf_module(periods
    src/periods/period_ranking.c src/periods/best_period_list.c)
target_link_libraries(trf_periods PUBLIC trf_core)

trf_module(alignment
    src/alignment/similarity.c src/alignment/newwrap.c
    src/alignment/narrowband.c src/alignment/narrowband_scalar.c
    src/alignment/narrowband_simd.c src/alignment/pair_alignment.c)
target_link_libraries(trf_alignment PUBLIC trf_core)

trf_module(output
    src/output/alignment_text.c src/output/html.c src/output/summary_html.c
    src/output/data_file.c src/output/mask.c)
target_link_libraries(trf_output PUBLIC trf_alignment trf_core)

trf_module(consensus
    src/consensus/consensus.c src/consensus/statistics.c)
target_link_libraries(trf_consensus PUBLIC trf_output trf_alignment trf_core)

trf_module(detection
    src/detection/tuple_scan.c src/detection/distance_list.c
    src/detection/distance_seen.c src/detection/criteria.c
    src/detection/coin_toss.c)
target_link_libraries(trf_detection PUBLIC
    trf_consensus trf_alignment trf_periods trf_core)

trf_module(redundancy
    src/redundancy/clean.c src/redundancy/list_sort.c src/redundancy/redundancy.c)
target_link_libraries(trf_redundancy PUBLIC trf_output trf_core)

trf_module(run
    src/run/control.c src/run/sequence.c src/run/fasta.c)
target_link_libraries(trf_run PUBLIC
    trf_detection trf_redundancy trf_output trf_periods trf_core)

add_executable(trf src/cli/main.c src/cli/arguments.c src/cli/usage.c)
target_link_libraries(trf PRIVATE trf_run trf_core)

# Install the program and, as the original TRF does, a link that carries the
# version and platform in its name.
include(GNUInstallDirs)
install(TARGETS trf RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(FILES README.md DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/doc/trf)
set(TRF_PLATFORM_NAME "")
if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|AMD64|amd64)$")
    set(TRF_WORD_SIZE 64)
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "^i.86$")
    set(TRF_WORD_SIZE 32)
else()
    set(TRF_WORD_SIZE "")
endif()
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    set(TRF_PLATFORM_NAME "linux${TRF_WORD_SIZE}")
elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    set(TRF_PLATFORM_NAME "mac")
elseif(CYGWIN OR MINGW)
    set(TRF_PLATFORM_NAME "dos${TRF_WORD_SIZE}")
endif()
if(UNIX)
    set(TRF_VERSIONED_NAME "trf${PROJECT_VERSION}.${TRF_PLATFORM_NAME}.exe")
    install(CODE "
        set(bindir \"\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}\")
        execute_process(COMMAND \${CMAKE_COMMAND} -E create_symlink
            trf${CMAKE_EXECUTABLE_SUFFIX} \"\${bindir}/${TRF_VERSIONED_NAME}\")
        message(STATUS \"Installing: \${bindir}/${TRF_VERSIONED_NAME}\")")
endif()

include(tests/CMakeLists.txt)
