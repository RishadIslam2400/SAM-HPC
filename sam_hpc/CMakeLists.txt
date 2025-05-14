find_package(Eigen3 REQUIRED)

add_library(sam INTERFACE)
target_include_directories(sam 
    INTERFACE 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../sam_hpc/include>
)
target_link_libraries(sam
    INTERFACE
    csr_matrix
    tools
    mgs_solver
    householder_solver
    Eigen3::Eigen
)