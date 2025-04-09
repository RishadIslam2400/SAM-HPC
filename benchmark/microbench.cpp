#include "cscMatrix.hpp"
#include "config.hpp"
#include "sparsityPattern.hpp"
#include "samFunc.hpp"
#include "samFunc_threads.hpp"

#include <iostream>
#include <vector>
#include <functional>
#include <chrono>
#include <string>
#include <fstream>
#include <sstream>

bool read_mat(const char *filename, csc_matrix<double, ptrdiff_t> &sparse_matrix)
{
    std::ifstream infile(filename);
    if (!infile.is_open())
    {
        std::cerr << "Error opening file: " << filename << std::endl;
        return false;
    }

    size_t nrows, ncols, nnz;
    std::string line;

    // First line: nrows ncols nnz
    if (!std::getline(infile, line))
        return false;
    std::istringstream header(line);
    if (!(header >> nrows >> ncols >> nnz))
    {
        std::cerr << "Error reading matrix dimensions" << std::endl;
        return false;
    }

    sparse_matrix.mNumRows = nrows;
    sparse_matrix.mNumCols = ncols;
    sparse_matrix.mNNZ = nnz;

    // Second line: column pointers
    if (!std::getline(infile, line))
        return false;
    std::istringstream colptr_stream(line);
    sparse_matrix.mColPointers.resize(ncols + 1);
    for (size_t i = 0; i <= ncols; ++i)
    {
        if (!(colptr_stream >> sparse_matrix.mColPointers[i]))
        {
            std::cerr << "Error reading column pointers" << std::endl;
            return false;
        }
    }

    // Third line: row indices
    if (!std::getline(infile, line))
        return false;
    std::istringstream rowind_stream(line);
    sparse_matrix.mRowIndices.resize(nnz);
    for (size_t i = 0; i < nnz; ++i)
    {
        if (!(rowind_stream >> sparse_matrix.mRowIndices[i]))
        {
            std::cerr << "Error reading row indices" << std::endl;
            return false;
        }
    }

    // Fourth line: values
    if (!std::getline(infile, line))
        return false;
    std::istringstream val_stream(line);
    sparse_matrix.mValues.resize(nnz);
    for (size_t i = 0; i < nnz; ++i)
    {
        if (!(val_stream >> sparse_matrix.mValues[i]))
        {
            std::cerr << "Error reading values" << std::endl;
            return false;
        }
    }

    infile.close();
    return true;
}

int main(int argc, char **argv)
{
    // get the configuration, print it
    config_t config;
    parseargs(argc, argv, config);
    config.print();

    // Create target matrix
    csc_matrix<> target_matrix;
    std::string target_filename = "/home/rishad/SAM-HPC/top_opt_matrices_small_csc/matrix_1.txt";
    if (!read_mat(target_filename.c_str(), target_matrix))
    {
        std::cerr << "Error reading the sparse matrix!" << std::endl;
        return 1;
    }

    // Read the reqired sparse matrices
    std::cout << "Sequential Benchmark: " << std::endl;
    std::cout << "SAM computation time (simple sparsity pattern):" << std::endl;
    for (int i = 1; i <= 50; i += 10)
    {
        std::string filename = "/home/rishad/SAM-HPC/top_opt_matrices_small_csc/matrix_" + std::to_string(i + 1) + ".txt";
        csc_matrix<> test_matrix;
        if (!read_mat(filename.c_str(), test_matrix))
        {
            std::cerr << "Error reading the sparse matrix!" << std::endl;
            return 1;
        }

        // Run the benchmark for simple sparsity pattern
        std::cout << "Matrix " << i << ": " << std::flush;
        std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
        for (int iter = 0; iter < config.iters; iter++)
        {
            sparsity_pattern<> sparsity_pattern;
            simple_sparsity_pattern(test_matrix, sparsity_pattern);

            // Compute the SAM
            csc_matrix<> map;
            SAM(test_matrix, target_matrix, sparsity_pattern, map);
        }
        // Compute average time
        std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        std::cout << elapsed.count() / config.iters << " s" << std::endl;
    }

    
    // test all other sparsity patterns for matrix_11
    std::string filename = "/home/rishad/SAM-HPC/top_opt_matrices_small_csc/matrix_11.txt";
    csc_matrix<> test_matrix;
    if (!read_mat(filename.c_str(), test_matrix))
    {
        std::cerr << "Error reading the sparse matrix!" << std::endl;
        return 1;
    }
    std::cout << "SAM computation time (global sparsity pattern):" << std::endl;
    
    // Run the benchmark for global sparsity pattern
    std::cout << "Matrix 11: " << std::flush;
    std::chrono::high_resolution_clock::time_point global_sparsity_pattern_start = std::chrono::high_resolution_clock::now();
    for (int iter = 0; iter < config.iters; iter++)
    {
        sparsity_pattern<> sparsity_pattern;
        sparsity_pattern_global_thresh(test_matrix, 0.001, sparsity_pattern);

        // Compute the SAM
        csc_matrix<> map;
        SAM(test_matrix, target_matrix, sparsity_pattern, map);
    }
    // Compute average time
    std::chrono::high_resolution_clock::time_point global_sparsity_pattern_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> global_sparsity_pattern_elapsed = global_sparsity_pattern_end - global_sparsity_pattern_start;
    std::cout << global_sparsity_pattern_elapsed.count() / config.iters << " s" << std::endl;

    // Run the benchmark for column sparsity pattern
    std::cout << "Matrix 11: " << std::flush;
    std::chrono::high_resolution_clock::time_point column_sparsity_pattern_start = std::chrono::high_resolution_clock::now();
    for (int iter = 0; iter < config.iters; iter++)
    {
        sparsity_pattern<> sparsity_pattern;
        sparsity_pattern_col_thresh(test_matrix, 0.8, sparsity_pattern);

        // Compute the SAM
        csc_matrix<> map;
        SAM(test_matrix, target_matrix, sparsity_pattern, map);
    }
    // Compute average time
    std::chrono::high_resolution_clock::time_point column_sparsity_pattern_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> column_sparsity_pattern_elapsed = column_sparsity_pattern_end - column_sparsity_pattern_start;
    std::cout << column_sparsity_pattern_elapsed.count() / config.iters << " s" << std::endl;

    // Run the benchmark for lfil sparsity pattern
    std::cout << "Matrix 11: " << std::flush;
    std::chrono::high_resolution_clock::time_point lfil_sparsity_pattern_start = std::chrono::high_resolution_clock::now();
    for (int iter = 0; iter < config.iters; iter++)
    {
        sparsity_pattern<> sparsity_pattern;
        sparsity_pattern_lfil_thresh(test_matrix, 5, sparsity_pattern);

        // Compute the SAM
        csc_matrix<> map;
        SAM(test_matrix, target_matrix, sparsity_pattern, map);
    }
    // Compute average time
    std::chrono::high_resolution_clock::time_point lfil_sparsity_pattern_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> lfil_sparsity_pattern_elapsed = lfil_sparsity_pattern_end - lfil_sparsity_pattern_start;
    std::cout << lfil_sparsity_pattern_elapsed.count() / config.iters << " s" << std::endl;

    // Parallel Benchmark
    // Read the reqired sparse matrices
    std::cout << std::endl << std::endl << "Parallel Benchmark: " << std::endl;
    std::cout << "SAM computation time (simple sparsity pattern):" << std::endl;
    for (int i = 1; i <= 50; i += 10)
    {
        std::string filename = "/home/rishad/SAM-HPC/top_opt_matrices_small_csc/matrix_" + std::to_string(i + 1) + ".txt";
        csc_matrix<> test_matrix;
        if (!read_mat(filename.c_str(), test_matrix))
        {
            std::cerr << "Error reading the sparse matrix!" << std::endl;
            return 1;
        }

        // Run the benchmark for simple sparsity pattern
        std::cout << "Matrix " << i << ": " << std::flush;
        std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
        for (int iter = 0; iter < config.iters; iter++)
        {
            sparsity_pattern<> sparsity_pattern;
            simple_sparsity_pattern(test_matrix, sparsity_pattern);

            // Compute the SAM
            csc_matrix<> map;
            SAM_std_thread(test_matrix, target_matrix, sparsity_pattern, map, config.threads);
        }
        // Compute average time
        std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        std::cout << elapsed.count() / config.iters << " s" << std::endl;
    }
    
    std::cout << "SAM computation time (global sparsity pattern):" << std::endl;

    // Run the benchmark for global sparsity pattern
    std::cout << "Matrix 11: " << std::flush;
    std::chrono::high_resolution_clock::time_point par_global_sparsity_pattern_start = std::chrono::high_resolution_clock::now();
    for (int iter = 0; iter < config.iters; iter++)
    {
        sparsity_pattern<> sparsity_pattern;
        sparsity_pattern_global_thresh(test_matrix, 0.001, sparsity_pattern);

        // Compute the SAM
        csc_matrix<> map;
        SAM_std_thread(test_matrix, target_matrix, sparsity_pattern, map, config.threads);
    }
    // Compute average time
    std::chrono::high_resolution_clock::time_point par_global_sparsity_pattern_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> par_global_sparsity_pattern_elapsed = par_global_sparsity_pattern_end - par_global_sparsity_pattern_start;
    std::cout << par_global_sparsity_pattern_elapsed.count() / config.iters << " s" << std::endl;

    // Run the benchmark for column sparsity pattern
    std::cout << "Matrix 11: " << std::flush;
    std::chrono::high_resolution_clock::time_point par_column_sparsity_pattern_start = std::chrono::high_resolution_clock::now();
    for (int iter = 0; iter < config.iters; iter++)
    {
        sparsity_pattern<> sparsity_pattern;
        sparsity_pattern_col_thresh(test_matrix, 0.8, sparsity_pattern);

        // Compute the SAM
        csc_matrix<> map;
        SAM_std_thread(test_matrix, target_matrix, sparsity_pattern, map, config.threads);
    }
    // Compute average time
    std::chrono::high_resolution_clock::time_point par_column_sparsity_pattern_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> par_column_sparsity_pattern_elapsed = par_column_sparsity_pattern_end - par_column_sparsity_pattern_start;
    std::cout << par_column_sparsity_pattern_elapsed.count() / config.iters << " s" << std::endl;

    // Run the benchmark for lfil sparsity pattern
    std::cout << "Matrix 11: " << std::flush;
    std::chrono::high_resolution_clock::time_point par_lfil_sparsity_pattern_start = std::chrono::high_resolution_clock::now();
    for (int iter = 0; iter < config.iters; iter++)
    {
        sparsity_pattern<> sparsity_pattern;
        sparsity_pattern_lfil_thresh(test_matrix, 5, sparsity_pattern);

        // Compute the SAM
        csc_matrix<> map;
        SAM_std_thread(test_matrix, target_matrix, sparsity_pattern, map, config.threads);
    }
    // Compute average time
    std::chrono::high_resolution_clock::time_point par_lfil_sparsity_pattern_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> par_lfil_sparsity_pattern_elapsed = par_lfil_sparsity_pattern_end - par_lfil_sparsity_pattern_start;
    std::cout << par_lfil_sparsity_pattern_elapsed.count() / config.iters << " s" << std::endl;

    return 0;
}