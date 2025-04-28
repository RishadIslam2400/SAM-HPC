cmake_minimum_required(VERSION 3.10)
project(SAM_HPC LANGUAGES CXX)

add_executable(qr_bench run_qr_bench.cpp)
target_link_libraries(qr_bench PUBLIC mgs_solver householder_solver sam)
target_compile_options(qr_bench PRIVATE -std=c++20 -Wall -Werror)