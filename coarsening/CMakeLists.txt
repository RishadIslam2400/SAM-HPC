add_library(coarsening INTERFACE)
target_include_directories(relaxation
    INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../coarsening/include>
)
target_link_libraries(coarsening INTERFACE csr_matrix solver)