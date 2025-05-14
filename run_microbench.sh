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

OUTPUT_FILE="top_opt_sam_result.txt"
SOURCE_FILE="top_opt_matrices_small_csr/matrix_6.txt"
TARGET_FILE="top_opt_matrices_small_csr/matrix_1.txt"
K=30

touch "$OUTPUT_FILE"
echo "Microbenchmark results" > "$OUTPUT_FILE"

parallel_benchmarks() {
    echo "----------Parallel----------" >> "$OUTPUT_FILE"
    for THREAD in 2 4 8 16 32 64 96 128 160
    do
        echo "Running with $THREAD threads..."
        echo "Threads: $THREAD" >> "$OUTPUT_FILE"
        ./build/benchmarks/microbench/microbench -n "sam microbenchmark parallel" -x "$SOURCE_FILE" -y "$TARGET_FILE" -k "$K" -t "$THREAD" >> "$OUTPUT_FILE"
    done

}

sequential_benchmark() {
    echo "----------Sequential----------" >> "$OUTPUT_FILE"
    echo "Running with 1 thread..."
    ./build/benchmarks/microbench/microbench -n "sam microbenchmark sequential" -x "$SOURCE_FILE" -y "$TARGET_FILE" -k "$K" -t 1 >> "$OUTPUT_FILE"
}

build_parallel
parallel_benchmarks

build_sequential
sequential_benchmark