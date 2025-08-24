add_library(preconditioner INTERFACE)
target_include_directories(preconditioner
    INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../preconditioner/include>
)
target_link_libraries(preconditioner INTERFACE relaxation coarsening)