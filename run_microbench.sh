#!/bin/bash

cd build
make clean
cmake -DCMAKE_BUILD_TYPE=Release ..
make
cd ..

# Output file
OUTPUT_FILE="microbench.txt"
echo "Microbenchmark Results" > "$OUTPUT_FILE"
echo "-----------Sequential-----------" >> "$OUTPUT_FILE"
./build/benchmark/microbench -n "microbenchmark" -f "/top_opt_matrices_small_csc/" -k 10 >> "$OUTPUT_FILE"

# Thread counts: geometric progression
echo "-----------Parallel-----------" >> "$OUTPUT_FILE"
for i in 2 4 8 32 64 92 128
do
    echo "Running with $i threads..."
    ./build/benchmark/microbenchmark_parallel -n "microbenchmark_parallel" -f "/top_opt_matrices_small_csc/" -k 10 -t "$i" >> "$OUTPUT_FILE"
done