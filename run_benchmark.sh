#!/bin/bash
set -e # Exit immediately if a command exits with a non-zero status.

# Configuration
# Use command-line arguments with defaults.
# Usage: ./run_benchmark.sh [MATRIX_DIR] [RHS_DIR] [K]
DEFAULT_MATRIX_DIR="/home/rishad/SAM-HPC/top_opt_matrices_small_csr/"
DEFAULT_RHS_DIR="/home/rishad/SAM-HPC/top_opt_rhs_small/"
DEFAULT_K=20

MATRIX_DIR="${1:-$DEFAULT_MATRIX_DIR}"
RHS_DIR="${2:-$DEFAULT_RHS_DIR}"
K="${3:-$DEFAULT_K}"

# Create a unique, timestamped output file to prevent overwriting results.
TIMESTAMP=$(date +%Y-%m-%d_%H-%M-%S)
OUTPUT_FILE="top_opt_results_${TIMESTAMP}.txt"
BENCHMARK_EXE="./build/benchmarks/top_opt_bench"
THREADS="1 2 4 8 16 32 64 96 128 160"

# Build the project only if the executable doesn't exist.
build_if_needed() {
    if [ ! -f "$BENCHMARK_EXE" ]; then
        echo "Building SAM HPC..."
        export CC=clang-20
        export CXX=clang++-20
        cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DENABLE_TESTS=off -DENABLE_BENCHMARKS=on
        cmake --build build
    else
        echo "Build artifacts found. Skipping build."
    fi
}

# Run the benchmarks for a given set of threads.
run_parallel_benchmarks() {
    printf "\n---------- Benchmarks ----------\n" >> "$OUTPUT_FILE"
    
    for THREAD in $THREADS; do
        echo "Running with $THREAD threads..."
        printf "\nThreads: %s\n" "$THREAD" >> "$OUTPUT_FILE"
        
        # Execute the benchmark and append its output.
        "$BENCHMARK_EXE" -n "sam with solver benchmark" -x "$MATRIX_DIR" -y "$RHS_DIR" -k "$K" -t "$THREAD" >> "$OUTPUT_FILE"
    done
}

# --- Main Execution ---
echo "Starting Topology Optimization Benchmark"
echo "Results will be saved to: $OUTPUT_FILE"

# Initialize the output file with a header.
{
    echo "Topology Optimization Benchmark Results"
    echo "Date: $(date)"
    echo "Matrix Directory: $MATRIX_DIR"
    echo "RHS Directory: $RHS_DIR"
    echo "K value: $K"
} > "$OUTPUT_FILE"


build_if_needed
run_parallel_benchmarks

echo "Benchmark finished successfully."
echo "Results available in $OUTPUT_FILE"