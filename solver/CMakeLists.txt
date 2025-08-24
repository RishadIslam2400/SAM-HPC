add_library(solver INTERFACE)
target_include_directories(solver
    INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../solver/include>
)
target_link_libraries(solver INTERFACE csr_matrix preconditioner)