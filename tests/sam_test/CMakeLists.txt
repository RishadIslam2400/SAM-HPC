cmake_minimum_required(VERSION 3.10)
project(SAM_HPC LANGUAGES CXX)

add_executable(
    sam_test
    run_sam_tests.cpp
)
target_link_libraries(sam_test PUBLIC csr_matrix sam csr_test_utils)
target_compile_options(sam_test PRIVATE -std=c++20 -Wall -Werror -fsanitize=address)
target_link_options(sam_test PRIVATE -fsanitize=address)