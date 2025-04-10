#!/bin/bash

cd build
make clean
cmake -DCMAKE_BUILD_TYPE=Release ..
make
cd ..

# Output file
OUTPUT_FILE="microbench.txt"
#echo "Microbenchmark Results" > "$OUTPUT_FILE"
#echo "-----------Sequential-----------" >> "$OUTPUT_FILE"
#./build/benchmark/microbench -n "microbenchmark" -f "/top_opt_matrices_small_csc/" -k 10 >> "$OUTPUT_FILE"

# Thread counts: geometric progression
echo "-----------Parallel-----------" >> "$OUTPUT_FILE"
echo >> "Running with 2 threads..."
./build/benchmark/microbench_parallel -n "microbenchmark_parallel" -f "/top_opt_matrices_small_csc/" -k 10 -t 2 >> "$OUTPUT_FILE"

echo >> "Running with 4 threads..."
./build/benchmark/microbench_parallel -n "microbenchmark_parallel" -f "/top_opt_matrices_small_csc/" -k 10 -t 4 >> "$OUTPUT_FILE"

echo >> "Running with 8 threads..."
./build/benchmark/microbench_parallel -n "microbenchmark_parallel" -f "/top_opt_matrices_small_csc/" -k 10 -t 8 >> "$OUTPUT_FILE"

echo >> "Running with 16 threads..."
./build/benchmark/microbench_parallel -n "microbenchmark_parallel" -f "/top_opt_matrices_small_csc/" -k 10 -t 16 >> "$OUTPUT_FILE"

echo >> "Running with 32 threads..."
./build/benchmark/microbench_parallel -n "microbenchmark_parallel" -f "/top_opt_matrices_small_csc/" -k 10 -t 32 >> "$OUTPUT_FILE"

echo >> "Running with 64 threads..."
./build/benchmark/microbench_parallel -n "microbenchmark_parallel" -f "/top_opt_matrices_small_csc/" -k 10 -t 64 >> "$OUTPUT_FILE"

echo >> "Running with 96 threads..."
./build/benchmark/microbench_parallel -n "microbenchmark_parallel" -f "/top_opt_matrices_small_csc/" -k 10 -t 96 >> "$OUTPUT_FILE"

echo >> "Running with 128 threads..."
./build/benchmark/microbench_parallel -n "microbenchmark_parallel" -f "/top_opt_matrices_small_csc/" -k 10 -t 128 >> "$OUTPUT_FILE"
