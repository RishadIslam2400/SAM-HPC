#include "config.hpp"
#include "sam.hpp"
#include "read_mat.hpp"
#include "amg.hpp"
#include "ruge_stuben.hpp"
#include "damped_jacobi.hpp"
#include "smoothed_aggregation.hpp"
#include "gauss_seidel.hpp"
#include "gmres.hpp"
#include "timer.hpp"

#include <map>

struct BenchmarkResults {
    uint64_t map_time_ms = 0;
    uint64_t solver_time_ms = 0;
    double relative_error = 0.0;
    size_t total_iterations = 0;
};

template <typename PatternType, typename... PatternArgs>
void run_benchmark(
    const config_t& cfg,
    const CSRMatrix<double>& targetMatrix,
    const amg<double, smoothed_aggregation, damped_jacobi>& P_0,
    const GMRES<double>::params& prm,
    std::map<std::string, BenchmarkResults>& results,
    const std::string& pattern_name,
    PatternArgs&&... args)
{
    std::cout << "Running benchmark: [" << pattern_name << "]..." << std::endl;

    for (int k = 0; k < cfg.iters; ++k) {
        for (int i = 1; i < 2; ++i) {
            // Data Loading (Not timed)
            std::string filepath = cfg.matrix_dir + "matrix_" + std::to_string(i + 1) + ".txt";
            CSRMatrix<double> A_i;
            read_mat<double>(filepath.c_str(), A_i);

            filepath = cfg.rhs_dir + "rhs_" + std::to_string(i + 1) + ".txt";
            std::vector<double> r_i = read_vec<double>(filepath.c_str());
            std::vector<double> x_i(A_i.m_cols, 0.0);

            CSRMatrix<double> map;
            Timer timer;

            // Time Map Computation
            timer.start();
            PatternType pattern_params(std::forward<PatternArgs>(args)...);
            SparsityPattern<double, PatternType> pattern(A_i, pattern_params);
            pattern.computePattern();
            SparseApproximateMap<double, PatternType>::computeMap(targetMatrix, A_i, pattern, map);
            results[pattern_name].map_time_ms += timer.elapsed();

            // Time Solver
            GMRES<double> solver(A_i.m_cols, prm);
            timer.start(); // Start timer right before the solve call
            auto [iters, error] = solver.solve(A_i, P_0, map, r_i, x_i);
            results[pattern_name].solver_time_ms += timer.elapsed();

            // Accumulate Stats
            results[pattern_name].total_iterations += iters;
            results[pattern_name].relative_error += error;
        }
    }
}

int main(int argc, char** argv) {
    config_t cfg;
    parseargs(argc, argv, cfg);
    cfg.print();

    num_threads = cfg.threads;
    tbb::global_control gc(tbb::global_control::max_allowed_parallelism, num_threads);

    // First matrix in the sequence, A_0 x_0 = r_0
    // read the matrix
    std::cout << "Matrix A" + std::to_string(0) + ':' << std::endl;
    std::string filepath = cfg.matrix_dir + "matrix_" + std::to_string(1) + ".txt";
    CSRMatrix<double> targetMatrix;
    read_mat<double>(filepath.c_str(), targetMatrix);

    // read the rhs vector
    filepath = cfg.rhs_dir + "rhs_" + std::to_string(1) + ".txt";
    std::vector<double> r_0 = read_vec<double>(filepath.c_str());
    std::vector<double> x_0(targetMatrix.m_cols, 0.0);

    // Compute the preconditioner
    amg<double, smoothed_aggregation, damped_jacobi> P_0(targetMatrix);
    GMRES<double>::params prm;
    prm.pside = precondSide::right;
    GMRES<double> solver_init(targetMatrix.m_cols, prm);

    // Solve the system
    auto [init_iters, init_error] = solver_init.solve(targetMatrix, P_0, r_0, x_0);
    std::cout << "Iterations: " << init_iters << ", Relative Error: " << init_error << "..." << std::endl;

    // There are 50 matrices in the sequence
    // Start from the second matrix and compute the map and solve the system
    // For now using the second matrix only

    // --- Benchmarking Section ---
    std::map<std::string, BenchmarkResults> results;

    // Call the generic benchmark function for each case.
    run_benchmark<SimplePattern>(cfg, targetMatrix, P_0, prm, results, "Simple Sparsity Pattern");

    run_benchmark<GlobalThresholdPattern>(cfg, targetMatrix, P_0, prm, results, "Global Sparsity (thresh=0.01)", 0.01);
    run_benchmark<GlobalThresholdPattern>(cfg, targetMatrix, P_0, prm, results, "Global Sparsity (thresh=0.001)", 0.001);
    run_benchmark<GlobalThresholdPattern>(cfg, targetMatrix, P_0, prm, results, "Global Sparsity (thresh=0.0001)", 0.0001);

    run_benchmark<ColumnThresholdPattern>(cfg, targetMatrix, P_0, prm, results, "Column Sparsity (thresh=0.7)", 0.7);
    run_benchmark<ColumnThresholdPattern>(cfg, targetMatrix, P_0, prm, results, "Column Sparsity (thresh=0.8)", 0.8);
    run_benchmark<ColumnThresholdPattern>(cfg, targetMatrix, P_0, prm, results, "Column Sparsity (thresh=0.9)", 0.9);

    run_benchmark<FixedNNZPattern>(cfg, targetMatrix, P_0, prm, results, "Fixed Sparsity (thresh = 3)", 3);
    run_benchmark<FixedNNZPattern>(cfg, targetMatrix, P_0, prm, results, "Fixed Sparsity (thresh = 5)", 5);
    run_benchmark<FixedNNZPattern>(cfg, targetMatrix, P_0, prm, results, "Fixed Sparsity (thresh = 7)", 7);


    // Report Final Results
    std::cout << "\n\nBenchmark Summary\n";
    std::cout << std::fixed << std::setprecision(4);

    for (const auto& pair : results) {
        const std::string& pattern_name = pair.first;
        const BenchmarkResults& res = pair.second;

        double avg_map_s = (res.map_time_ms / cfg.iters) / 1000000.0;
        double avg_solver_s = (res.solver_time_ms / cfg.iters) / 1000000.0;
        double avg_total_s = avg_map_s + avg_solver_s;
        double avg_iters = res.total_iterations / cfg.iters;
        double avg_relative_error = res.relative_error / cfg.iters;

        std::cout << "\nPattern: [" << pattern_name << "]\n";
        std::cout << "  - Average Map Time:    " << avg_map_s << " ms\n";
        std::cout << "  - Average Solver Time: " << avg_solver_s << " ms\n";
        std::cout << "  - Average Total Time:  " << avg_total_s << " ms\n";
        std::cout << "  - Average Iterations:  " << avg_iters << "\n";
        std::cout << "  - Average Rel Error:   " << std::scientific << avg_relative_error << std::fixed << "\n";
    }

}