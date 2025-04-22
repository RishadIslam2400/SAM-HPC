cmake_minimum_required(VERSION 3.10)
project(SAM_HPC LANGUAGES CXX)

add_executable(
    csr_test
    run_csr_tests.cpp   
)
target_link_libraries(csr_test PUBLIC csr_test_utils)
target_compile_options(csr_test PRIVATE -std=c++20 -Wall -Werror -fsanitize=address)
target_link_options(csr_test PRIVATE -fsanitize=address)