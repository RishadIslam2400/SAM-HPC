project(${CMAKE_PROJECT_NAME})
cmake_minimum_required(VERSION 3.10)

# set(MATLAB_ROOT "/usr/local/MATLAB/R2022b")
# include_directories(${MATLAB_ROOT}/extern/include)
# link_directories(${MATLAB_ROOT}/bin/glnxa64)

add_library(csc_matrix INTERFACE cscMatrix.hpp spgemm.hpp)
target_include_directories(csc_matrix INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/../include/)

add_library(sam_func INTERFACE samFunc.hpp sparsityPattern.hpp cscMatrix.hpp spgemm.hpp)
target_include_directories(sam_func INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/../include/)
target_link_libraries(sam_func INTERFACE householderQR mgsQR)

add_subdirectory(householder)
add_subdirectory(mgs)
add_subdirectory(config)