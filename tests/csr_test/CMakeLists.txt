cmake_minimum_required(VERSION 3.10)
project(SAM_HPC LANGUAGES CXX)

add_executable(
    csr_test
    run_csr_tests.cpp   
)
target_link_libraries(csr_test PUBLIC csr_test_utils)