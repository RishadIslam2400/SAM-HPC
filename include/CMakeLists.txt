project(${CMAKE_PROJECT_NAME})
cmake_minimum_required(VERSION 3.10)

set(MATLAB_ROOT "/usr/local/MATLAB/R2022b")
include_directories(${MATLAB_ROOT}/extern/include)
link_directories(${MATLAB_ROOT}/bin/glnxa64)

add_library(sam_func INTERFACE SAM_Func.hpp)
target_include_directories(sam_func INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/../include)
target_link_libraries(sam_func INTERFACE householderQR)

add_subdirectory(householder)
add_subdirectory(mgs)