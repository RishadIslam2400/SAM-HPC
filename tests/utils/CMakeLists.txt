add_library(csr_test_utils INTERFACE)
target_include_directories(csr_test_utils
    INTERFACE 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../utils/>
)
target_link_libraries(csr_test_utils INTERFACE csr_matrix)