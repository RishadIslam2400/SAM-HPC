#!/bin/bash
set -e

build_parallel() {
    echo "Building SAM Parallel..."
    rm -rf build
    mkdir -p build
    cd build

    cmake -DCMAKE_BUILD_TYPE=Release -DENABLE_TESTS=off -DENABLE_BENCHMARKS=on -DENABLE_SEQUENTIAL=off ..
    make
    cd .. 
}

build_sequential() {
    echo "Building SAM Sequential..."
    rm -rf build
    mkdir -p build
    cd build

    cmake -DCMAKE_BUILD_TYPE=Release -DENABLE_TESTS=off -DENABLE_BENCHMARKS=on -DENABLE_SEQUENTIAL=on ..
    make
    cd .. 
}

OUTPUT_FILE="demo_result.txt"
MATRIX_DIR="/home/rishad/SAM-HPC/top_opt_matrices_small_csr/"
RHS_DIR="/home/rishad/SAM-HPC/top_opt_matrices_small_csr/"
K=20

touch "$OUTPUT_FILE"
echo "demo results" > "$OUTPUT_FILE"

parallel_benchmarks() {
    echo "----------Parallel----------" >> "$OUTPUT_FILE"
    for THREAD in 2 4 8 16 32 64 96 128 160
    do
        echo "Running with $THREAD threads..."
        echo "Threads: $THREAD" >> "$OUTPUT_FILE"
        ./build/benchmarks/microbench/microbench -n "column threshold benchmark" -x "$MATRIX_DIR" -y "$RHS_DIR" -k "$K" -t "$THREAD" >> "$OUTPUT_FILE"
    done

}

build_parallel
parallel_benchmarks