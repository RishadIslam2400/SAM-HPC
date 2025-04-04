cmake_minimum_required(VERSION 3.10)
project(SAM_HPC LANGUAGES CXX)

add_subdirectory(utils)

add_executable(
    csr_test
    run_tests.cpp   
)
target_link_libraries(csr_test PUBLIC csr_test_utils)