add_library(householder_solver INTERFACE)
target_include_directories(householder_solver
    INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../householder_solver/include>
)