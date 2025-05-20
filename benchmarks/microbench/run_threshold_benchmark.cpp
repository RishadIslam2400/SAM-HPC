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

    CSRMatrix<double> source = read_mat<double>(cfg.source_file.c_str());
    std::cout << "Source matrix: A5" << std::endl;

    // Benchmarking SAM with Column threshold with different threhold parameters
    // This is comparing the total SAM computation
    Timer timer;
    std::cout << "Column threshold parameter - 0.8 (SAM computation): " << std::flush;
    for (int i = 0; i < cfg.iters; ++i) {
        // Column/Row sparsity pattern
        ColumnThresholdPattern thresh{0.8}; // arbitrary
        SparsityPattern<double, ColumnThresholdPattern> pattern(source, thresh);
        pattern.computePattern();

        // Compute map
        SparseApproximateMap<double, ColumnThresholdPattern> map(target, source, pattern);
        map.computeMap();
    }
    std::cout << (timer.elapsed() / cfg.iters) / 1000000.0 << " s" << std::endl;

    std::cout << "Column/Row threhold parameter - 0.8 (Sparsity pattern computation): " << std::flush;
    timer.start();
    for (int i = 0; i < cfg.iters; ++i) {
        // Column/Row sparsity pattern
        ColumnThresholdPattern thresh{0.8}; // arbitrary
        SparsityPattern<double, ColumnThresholdPattern> pattern(source, thresh);
        pattern.computePattern();
    }
    std::cout << (timer.elapsed() / cfg.iters) / 1000000.0 << " s" << std::endl;

    std::cout << "Column/Row Sparsity Pattern (nnz): " << std::flush;
    ColumnThresholdPattern col_thresh1{0.8};
    SparsityPattern<double, ColumnThresholdPattern> columnPattern1(source, col_thresh1);
    columnPattern1.computePattern();
    std::cout << columnPattern1.getNNZ() << std::endl;

    std::cout << "Column threshold parameter - 0.85 (SAM computation): " << std::flush;
    for (int i = 0; i < cfg.iters; ++i) {
        // Column/Row sparsity pattern
        ColumnThresholdPattern thresh{0.85}; // arbitrary
        SparsityPattern<double, ColumnThresholdPattern> pattern(source, thresh);
        pattern.computePattern();

        // Compute map
        SparseApproximateMap<double, ColumnThresholdPattern> map(target, source, pattern);
        map.computeMap();
    }
    std::cout << (timer.elapsed() / cfg.iters) / 1000000.0 << " s" << std::endl;

    std::cout << "Column/Row threhold parameter - 0.85 (Sparsity pattern computation): " << std::flush;
    timer.start();
    for (int i = 0; i < cfg.iters; ++i) {
        // Column/Row sparsity pattern
        ColumnThresholdPattern thresh{0.85}; // arbitrary
        SparsityPattern<double, ColumnThresholdPattern> pattern(source, thresh);
        pattern.computePattern();
    }
    std::cout << (timer.elapsed() / cfg.iters) / 1000000.0 << " s" << std::endl;

    std::cout << "Column/Row Sparsity Pattern (nnz): " << std::flush;
    ColumnThresholdPattern col_thresh2{0.85};
    SparsityPattern<double, ColumnThresholdPattern> columnPattern2(source, col_thresh2);
    columnPattern2.computePattern();
    std::cout << columnPattern2.getNNZ() << std::endl;

    std::cout << "Column threshold parameter - 0.9 (SAM computation): " << std::flush;
    for (int i = 0; i < cfg.iters; ++i) {
        // Column/Row sparsity pattern
        ColumnThresholdPattern thresh{0.9}; // arbitrary
        SparsityPattern<double, ColumnThresholdPattern> pattern(source, thresh);
        pattern.computePattern();

        // Compute map
        SparseApproximateMap<double, ColumnThresholdPattern> map(target, source, pattern);
        map.computeMap();
    }
    std::cout << (timer.elapsed() / cfg.iters) / 1000000.0 << " s" << std::endl;

    std::cout << "Column/Row threhold parameter - 0.9 (Sparsity pattern computation): " << std::flush;
    timer.start();
    for (int i = 0; i < cfg.iters; ++i) {
        // Column/Row sparsity pattern
        ColumnThresholdPattern thresh{0.9}; // arbitrary
        SparsityPattern<double, ColumnThresholdPattern> pattern(source, thresh);
        pattern.computePattern();
    }
    std::cout << (timer.elapsed() / cfg.iters) / 1000000.0 << " s" << std::endl;

    std::cout << "Column/Row Sparsity Pattern (nnz): " << std::flush;
    ColumnThresholdPattern col_thresh3{0.9};
    SparsityPattern<double, ColumnThresholdPattern> columnPattern3(source, col_thresh3);
    columnPattern3.computePattern();
    std::cout << columnPattern3.getNNZ() << std::endl;

    std::cout << "Column threshold parameter - 0.95 (SAM computation): " << std::flush;
    for (int i = 0; i < cfg.iters; ++i) {
        // Column/Row sparsity pattern
        ColumnThresholdPattern thresh{0.95}; // arbitrary
        SparsityPattern<double, ColumnThresholdPattern> pattern(source, thresh);
        pattern.computePattern();

        // Compute map
        SparseApproximateMap<double, ColumnThresholdPattern> map(target, source, pattern);
        map.computeMap();
    }
    std::cout << (timer.elapsed() / cfg.iters) / 1000000.0 << " s" << std::endl;

    std::cout << "Column/Row threhold parameter - 0.95 (Sparsity pattern computation): " << std::flush;
    timer.start();
    for (int i = 0; i < cfg.iters; ++i) {
        // Column/Row sparsity pattern
        ColumnThresholdPattern thresh{0.95}; // arbitrary
        SparsityPattern<double, ColumnThresholdPattern> pattern(source, thresh);
        pattern.computePattern();
    }
    std::cout << (timer.elapsed() / cfg.iters) / 1000000.0 << " s" << std::endl;

    std::cout << "Column/Row Sparsity Pattern (nnz): " << std::flush;
    ColumnThresholdPattern col_thresh4{0.95};
    SparsityPattern<double, ColumnThresholdPattern> columnPattern4(source, col_thresh4);
    columnPattern4.computePattern();
    std::cout << columnPattern4.getNNZ() << std::endl;

    std::cout << "Column threshold parameter - 0.99 (SAM computation): " << std::flush;
    for (int i = 0; i < cfg.iters; ++i) {
        // Column/Row sparsity pattern
        ColumnThresholdPattern thresh{0.99}; // arbitrary
        SparsityPattern<double, ColumnThresholdPattern> pattern(source, thresh);
        pattern.computePattern();

        // Compute map
        SparseApproximateMap<double, ColumnThresholdPattern> map(target, source, pattern);
        map.computeMap();
    }
    std::cout << (timer.elapsed() / cfg.iters) / 1000000.0 << " s" << std::endl;

    std::cout << "Column/Row threhold parameter - 0.99 (Sparsity pattern computation): " << std::flush;
    timer.start();
    for (int i = 0; i < cfg.iters; ++i) {
        // Column/Row sparsity pattern
        ColumnThresholdPattern thresh{0.99}; // arbitrary
        SparsityPattern<double, ColumnThresholdPattern> pattern(source, thresh);
        pattern.computePattern();
    }
    std::cout << (timer.elapsed() / cfg.iters) / 1000000.0 << " s" << std::endl;

    std::cout << "Column/Row Sparsity Pattern (nnz): " << std::flush;
    ColumnThresholdPattern col_thresh5{0.99};
    SparsityPattern<double, ColumnThresholdPattern> columnPattern5(source, col_thresh5);
    columnPattern5.computePattern();
    std::cout << columnPattern5.getNNZ() << std::endl;
}