#!/bin/sh
# Build the original TRF for the differential tests:
#   <output>/source.tar   original sources (TRF 4.10.0, commit 355c1f9)
#   <output>/trf          original executable, built with its default flags
# Usage: tools/oracle/build_original_trf.sh [output]   (default: build/oracle)
# Needs a clone that includes the upstream history, a C compiler, and make.
# Environment: JOBS (default 4), CC, CFLAGS.
set -eu
commit=355c1f92b24f6b2e71650eda3549d46f0316e571
root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
output=${1:-"$root/build/oracle"}
mkdir -p "$output"
output=$(CDPATH= cd -- "$output" && pwd)
git -C "$root" archive --format=tar -o "$output/source.tar" "$commit"
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
tar -xf "$output/source.tar" -C "$work"
(
    cd "$work"
    # The original ships generated Autotools files. Touching them in
    # dependency order keeps make from trying to regenerate them.
    touch aclocal.m4
    touch configure
    touch config.h.in
    touch Makefile.in src/Makefile.in
    mkdir build
    cd build
    ../configure > configure.log 2>&1
    make -j "${JOBS:-4}" > make.log 2>&1
)
cp "$work/build/src/trf" "$output/trf"
echo "Prepared: $output/source.tar and $output/trf"
