#include "cases/test_damped_jacobi.cpp"
#include "cases/test_gauss_seidel.cpp"
#include "cases/test_aggregation.cpp"
#include "cases/test_ruge_stuben.cpp"
#include "cases/test_cuthill_mckee.cpp"
#include "cases/test_skyline_lu.cpp"
#include "cases/test_amg.cpp"
#include "cases/test_gmres.cpp"

int main() {
    std::cout << "Running Damped Jacobi Relaxation test..." << std::endl;
    test_jacobi_constructor();
    test_jacobi_presmoothing();
    test_jacobi_postsmoothing();
    test_direct_jacobi();

    std::cout << "\nRunning Gauss Seidel Relaxation test..." << std::endl;
    test_gauss_siedel_presmoothing_serial();
    test_gauss_siedel_postsmoothing_serial();
    test_gauss_siedel_direct_sovler_serial();
    test_gauss_siedel_presmoothing_parallel();
    test_gauss_siedel_postsmoothing_parallel();
    test_gauss_siedel_direct_sovler_parallel();

    std::cout << "\nRunning Smoothed Aggregation Coarsening test..." << std::endl;
    test_basic_aggregation();
    test_multiple_aggregates();
    test_remove_small_aggregates();
    test_piecewise_constant_prolongation();
    test_nullspace_prolongation();
    test_smoothed_transfer_operator();

    std::cout << "\nRunning Ruge-Stuben C/F Splitting test..." << std::endl;
    test_1d_laplacian_splitting_and_interpolation();

    std::cout << "\nRunning AMG test..." << std::endl;
    test_simple_path_graph_cm();
    test_reverse_ordering_rcm();
    test_disconnected_graph();
    test_skyline_lu_solver();
    test_standard_v_cycle_1();
    test_standard_v_cycle_2();

    std::cout << "\nRunning GMRES solver test..." << std::endl;
    test_gmres_no_precond();
    test_gmres_amg_precond_1();
    test_gmres_amg_precond_2();
}