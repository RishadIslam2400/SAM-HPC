cmake_minimum_required(VERSION 3.10)
project(SAM_HPC LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)

set(MATLAB_ROOT "/usr/local/MATLAB/R2022b")
include_directories(${MATLAB_ROOT}/extern/include)
link_directories(${MATLAB_ROOT}/bin/glnxa64)

add_subdirectory(include)

add_executable(test_matio test_matio.cpp)
target_link_libraries(test_matio PUBLIC mat mx)

# option(BUILD_EXAMPLES "Build examples" ON)
# if(BUILD_EXAMPLES)
# message(STATUS "Building examples")
# add_subdirectory(examples)
# endif()
add_subdirectory(examples)