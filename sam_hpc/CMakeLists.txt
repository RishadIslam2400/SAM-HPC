project(${CMAKE_PROJECT_NAME})
cmake_minimum_required(VERSION 3.10)

add_library(sam INTERFACE)
target_include_directories(sam INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/../sam_hpc/include)
target_link_libraries(sam INTERFACE csr_matrix)
target_compile_features(sam INTERFACE cxx_std_20)

add_subdirectory(src)