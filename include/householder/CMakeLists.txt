project(${CMAKE_PROJECT_NAME})

find_package(TBB REQUIRED COMPONENTS tbb)
add_library(householderQR INTERFACE householderQR.hpp)
target_include_directories(householderQR INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/../householder)
target_link_libraries(householderQR INTERFACE tbb)