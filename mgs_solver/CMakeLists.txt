add_library(mgs_solver INTERFACE)
target_include_directories(mgs_solver
    INTERFACE 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../mgs_solver/include/>
)