#!/bin/bash

cd build
make clean
cmake ..
make
cd ..

# Output file
OUTPUT_FILE="microbench.txt"
echo "Microbenchmark Results" > "$OUTPUT_FILE"

# Thread counts: geometric progression
for i in 2 4 8 32 64 92 128
do
    echo "Running with $i threads..."
    ./build/benchmark/microbench -n "microbenchmark" -f "/top_opt_matrices_small_csc/" -k 10 -t "$i" >> "$OUTPUT_FILE"
done