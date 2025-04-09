cmake_minimum_required(VERSION 3.10)
project(SAM_HPC LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)

# set(MATLAB_ROOT "/usr/local/MATLAB/R2022b")
# include_directories(${MATLAB_ROOT}/extern/include)
# link_directories(${MATLAB_ROOT}/bin/glnxa64)

add_subdirectory(include)

# option(BUILD_EXAMPLES "Build examples" ON)
# if(BUILD_EXAMPLES)
# message(STATUS "Building examples")
# add_subdirectory(examples)
# endif()
add_subdirectory(examples)
add_subdirectory(benchmark)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")