project(${CMAKE_PROJECT_NAME})
cmake_minimum_required(VERSION 3.10)

add_library(sam_func INTERFACE SAM_Func.hpp)
target_include_directories(sam_func INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/../include)