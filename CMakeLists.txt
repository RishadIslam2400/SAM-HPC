project(SAM_HPC)
cmake_minimum_required(VERSION 3.10)

set(MATLAB_ROOT "/usr/local/MATLAB/R2022b")
include_directories(${MATLAB_ROOT}/extern/include)
link_directories(${MATLAB_ROOT}/bin/glnxa64)

add_executable(test_matio test_matio.cpp)
target_link_libraries(test_matio PUBLIC mat mx)
