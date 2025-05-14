add_library(csr_matrix INTERFACE)
target_include_directories(csr_matrix
    INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../csr_matrix/include>
)

find_package(TBB REQUIRED)
find_package(Threads REQUIRED)
target_link_libraries(csr_matrix INTERFACE TBB::tbb Threads::Threads tools)