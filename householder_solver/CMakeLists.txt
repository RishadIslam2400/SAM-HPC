project(SAM_HPC LANGUAGES CXX)
cmake_minimum_required(VERSION 3.10)

add_library(householder_solver INTERFACE)
target_include_directories(householder_solver INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/../householder_solver/include/)