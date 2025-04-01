cmake_minimum_required(VERSION 3.10)
project(SAM_HPC LANGUAGES CXX)
set (CMAKE_CXX_STANDARD 20)

add_subdirectory(sam_hpc)
add_subdirectory(csr_matrix)
add_subdirectory(tests)