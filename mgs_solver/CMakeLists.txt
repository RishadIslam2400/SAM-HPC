project(SAM_HPC LANGUAGES CXX)
cmake_minimum_required(VERSION 3.10)

add_library(mgs_solver INTERFACE)
target_include_directories(mgs_solver INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/../mgs_solver/include/)