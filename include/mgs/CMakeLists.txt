project(${CMAKE_PROJECT_NAME})

find_package(TBB REQUIRED COMPONENTS tbb)
add_library(mgsQR INTERFACE mgsQR.hpp)
target_include_directories(mgsQR INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/../mgs)
target_link_libraries(mgsQR INTERFACE tbb)