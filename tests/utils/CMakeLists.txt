add_library(test_utils INTERFACE)
target_include_directories(test_utils
    INTERFACE 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../utils/>
)
target_link_libraries(test_utils INTERFACE csr_matrix)