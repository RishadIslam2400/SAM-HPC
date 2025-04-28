cmake_minimum_required(VERSION 3.10)
project(SAM_HPC LANGUAGES CXX)

add_executable(
    qr_test
    run_qr_tests.cpp
)
target_link_libraries(qr_test PUBLIC mgs_solver householder_solver sam csr_test_utils)
target_compile_options(qr_test PRIVATE -std=c++20 -Wall -Werror)