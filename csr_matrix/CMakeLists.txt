project(SAM_HPC LANGUAGES CXX)
cmake_minimum_required(VERSION 3.10)

add_library(csr_matrix INTERFACE)
target_include_directories(csr_matrix INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/../csr_matrix/include/)
target_compile_features(csr_matrix INTERFACE cxx_std_20)