#include "sam.hpp"
#include "config.hpp"
#include "read_mat.hpp"
#include "timer.hpp"

int main(int argc, char** argv) {
    config_t cfg;
    parseargs(argc, argv, cfg);
    cfg.print();

    num_threads = cfg.threads;
    
    // Parallel benchmarking
    // Mapping the fifth matrix A5 to the initial matrix A0
    // Read the target matrix
    CSRMatrix<double> target = read_mat<double>(cfg.target_file.c_str());
    std::cout << "Target matrix : A0" << std::endl;

    // std::string source_path = "../../../top_opt_small_csr/matrix_6.txt";
    // this should be the relative path of the input file to be included in the bash script
    CSRMatrix<double> source = read_mat<double>(cfg.source_file.c_str());
    std::cout << "Source matrix: A5" << std::endl;

    // Benchmarking SAM for 4 different sparsity patterns
    // This is comparing the total SAM computation
    std::cout << "Simple Sparsity Pattern (SAM computation): " << std::flush;
    Timer simple_timer{};
    simple_timer.start();
    for (int i = 0; i < cfg.iters; ++i) {
        // Compute sparsity pattern
        SparsityPattern<double, SimplePattern> pattern(source, SimplePattern{});
        pattern.computePattern();

        // Compute map
        SparseApproximateMap<double, SimplePattern> map(target, source, pattern);
        map.computeMap();
    }
    std::cout << (simple_timer.elapsed() / cfg.iters) / 1000000.0 << " s" << std::endl;

    std::cout << "Global Sparsity Pattern (SAM computation): " << std::flush;
    Timer global_timer{};
    global_timer.start();
    for (int i = 0; i < cfg.iters; ++i) {
        // Global sparsity pattern
        GlobalThresholdPattern thresh{0.001}; // arbitrary
        SparsityPattern<double, GlobalThresholdPattern> pattern(source, thresh);
        pattern.computePattern();

        // Compute map
        SparseApproximateMap<double, GlobalThresholdPattern> map(target, source, pattern);
        map.computeMap();
    }
    std::cout << (global_timer.elapsed() / cfg.iters) / 1000000.0 << " s" << std::endl;

    std::cout << "Column/Row Sparsity Pattern (SAM computation): " << std::flush;
    Timer local_timer{};
    local_timer.start();
    for (int i = 0; i < cfg.iters; ++i) {
        // Column/Row sparsity pattern
        ColumnThresholdPattern thresh{0.8}; // arbitrary
        SparsityPattern<double, ColumnThresholdPattern> pattern(source, thresh);
        pattern.computePattern();

        // Compute map
        SparseApproximateMap<double, ColumnThresholdPattern> map(target, source, pattern);
        map.computeMap();
    }
    std::cout << (local_timer.elapsed() / cfg.iters) / 1000000.0 << " s" << std::endl;

    std::cout << "Fixed Sparsity Pattern (SAM computation): " << std::flush;
    Timer fixed_timer{};
    fixed_timer.start();
    for (int i = 0; i < cfg.iters; ++i) {
        // Fixed NNZ sparsity pattern
        FixedNNZPattern thresh{5}; // arbitrary
        SparsityPattern<double, FixedNNZPattern> pattern(source, thresh);
        pattern.computePattern();

        // Compute map
        SparseApproximateMap<double, FixedNNZPattern> map(target, source, pattern);
        map.computeMap();
    }
    std::cout << (fixed_timer.elapsed() / cfg.iters) / 1000000.0 << " s" << std::endl;

    // This is comparing the sparsity pattern computation
    std::cout << "Simple Sparsity Pattern (Sparsity pattern computation): " << std::flush;
    simple_timer.start();
    for (int i = 0; i < cfg.iters; ++i) {
        // Compute sparsity pattern
        SparsityPattern<double, SimplePattern> pattern(source, SimplePattern{});
        pattern.computePattern();
    }
    std::cout << (simple_timer.elapsed() / cfg.iters) / 1000000.0 << " s" << std::endl;

    std::cout << "Global Sparsity Pattern (Sparsity pattern computation): " << std::flush;
    global_timer.start();
    for (int i = 0; i < cfg.iters; ++i) {
        // Global sparsity pattern
        GlobalThresholdPattern thresh{0.001}; // arbitrary
        SparsityPattern<double, GlobalThresholdPattern> pattern(source, thresh);
        pattern.computePattern();
    }
    std::cout << (global_timer.elapsed() / cfg.iters) / 1000000.0 << " s" << std::endl;

    std::cout << "Column/Row Sparsity Pattern (Sparsity pattern computation): " << std::flush;
    local_timer.start();
    for (int i = 0; i < cfg.iters; ++i) {
        // Column/Row sparsity pattern
        ColumnThresholdPattern thresh{0.8}; // arbitrary
        SparsityPattern<double, ColumnThresholdPattern> pattern(source, thresh);
        pattern.computePattern();
    }
    std::cout << (local_timer.elapsed() / cfg.iters) / 1000000.0 << " s" << std::endl;

    std::cout << "Fixed Sparsity Pattern (Sparsity pattern computation): " << std::flush;
    fixed_timer.start();
    for (int i = 0; i < cfg.iters; ++i) {
        // Fixed NNZ sparsity pattern
        FixedNNZPattern thresh{5}; // arbitrary
        SparsityPattern<double, FixedNNZPattern> pattern(source, thresh);
        pattern.computePattern();
    }
    std::cout << (fixed_timer.elapsed() / cfg.iters) / 1000000.0 << " s" << std::endl;

    // Count the number of non zero for each pattern
    std::cout << "Simple Sparsity Pattern (nnz): " << std::flush;
    SparsityPattern<double, SimplePattern> simplePattern(source, SimplePattern{});
    simplePattern.computePattern();
    std::cout << simplePattern.getNNZ() << std::endl;

    std::cout << "Global Sparsity Pattern (nnz): " << std::flush;
    GlobalThresholdPattern global_thresh{0.001};
    SparsityPattern<double, GlobalThresholdPattern> globalPattern(source, global_thresh);
    globalPattern.computePattern();
    std::cout << globalPattern.getNNZ() << std::endl;

    std::cout << "Column/Row Sparsity Pattern (nnz): " << std::flush;
    ColumnThresholdPattern col_thresh{0.8};
    SparsityPattern<double, ColumnThresholdPattern> columnPattern(source, col_thresh);
    columnPattern.computePattern();
    std::cout << columnPattern.getNNZ() << std::endl;

    std::cout << "Fixed Sparsity Pattern (nnz): " << std::flush;
    FixedNNZPattern lfil{5};
    SparsityPattern<double, FixedNNZPattern> fixedPattern(source, lfil);
    fixedPattern.computePattern();
    std::cout << fixedPattern.getNNZ() << std::endl;
}