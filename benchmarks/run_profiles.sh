#!/usr/bin/env bash
# SPDX-FileCopyrightText: 2026 Albert M. W. Yau
# SPDX-License-Identifier: GPL-2.0-or-later

# Reproduce the three dispatch profiles on the same runner.
set -euo pipefail

output_root="${1:?usage: run_profiles.sh OUTPUT_DIRECTORY}"
mkdir -p "$output_root"
conda list -n pyspharm-reference --explicit > "$output_root/conda-explicit.txt"

# Reject unsupported hosts rather than silently benchmarking a lower SIMD tier.
for feature in avx512f avx512dq avx512bw avx512vl; do
  grep -qm1 -w "$feature" /proc/cpuinfo || {
    echo "Missing $feature; an x86-64-v4 benchmark requires an AVX-512 runner" >&2
    exit 1
  }
done

for cap in 1 3 4; do
  export DUCC0_MAX_PSABI_LEVEL="$cap"
  output="$output_root/v$cap"
  mkdir -p "$output"

  conda run --no-capture-output -n pyspharm-reference python benchmarks/benchmark.py \
    --expected-profile "x86-64-v$cap" --threads 0 --warmups 2 --repeats 7 \
    --output "$output/multi-thread.csv"
  conda run --no-capture-output -n pyspharm-reference python benchmarks/benchmark.py \
    --expected-profile "x86-64-v$cap" --cases cc-2p5-1t \
    --threads 1 --warmups 2 --repeats 7 \
    --output "$output/one-thread.csv"
  conda run --no-capture-output -n pyspharm-reference python benchmarks/plot.py \
    --csv "$output/multi-thread.csv" "$output/one-thread.csv" \
    --output-dir "$output/figures"
done
