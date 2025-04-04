cmake_minimum_required(VERSION 3.10)
project(SAM_HPC LANGUAGES CXX)

add_library(csr_test_utils INTERFACE)
target_include_directories(csr_test_utils INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/../utils/)
target_link_libraries(csr_test_utils INTERFACE csr_matrix)
target_compile_features(csr_test_utils INTERFACE cxx_std_20)