#include "cscMatrix.hpp"
#include "config.hpp"
#include "sparsityPattern.hpp"
#include "samFunc.hpp"

#include <mat.h>
#include <matrix.h>
#include <iostream>
#include <vector>
#include <functional>
#include <chrono>

// This benchmark will compare the performance between all the different sparsity pattern choices
// We will compare the different sparsity pattern choices as a whole
// and then look into each spasity pattern to see which part of the code is taking the most time

bool read_mat(const char *filename, csc_matrix<double, ptrdiff_t> &sparse_matrix)
{
    // Read the .mat file
    MATFile *matfile = matOpen(filename, "r");
    if (!matfile)
    {
        std::cerr << "Error opening file: " << filename << std::endl;
        return false;
    }

    // Get all variable names in the file
    int num_variables;
    const char **dir = (const char **)matGetDir(matfile, &num_variables);
    if (!dir)
    {
        std::cerr << "Error getting variables" << std::endl;
        matClose(matfile);
        return false;
    }

    // Read the sparse matrix variable
    mxArray *sparseArray = nullptr;

    // Loop through the variables to find a sparse matrix
    for (int i = 0; i < num_variables; ++i)
    {
        mxArray *var = matGetVariable(matfile, dir[i]);
        if (!var)
        {
            std::cerr << "Error getting variable: " << dir[i] << std::endl;
            mxFree(dir);
            matClose(matfile);
            return false;
        }

        // Check if the variable is a sparse matrix
        if (mxIsSparse(var))
        {
            sparseArray = var; // Found the sparse matrix
            break;             // Stop searching
        }

        mxDestroyArray(var); // Free memory for non-sparse variables
    }
    mxFree(dir);

    // Ensure the variable is a sparse matrix
    if (!sparseArray)
    {
        std::cerr << "No sparse matrix found in the file!" << std::endl;
        matClose(matfile);
        return false;
    }

    // Get the data - Matlab stores sparse matrix in csc format
    double *val = mxGetPr(sparseArray);          // Non-zero elements
    mwIndex *rowIndices = mxGetIr(sparseArray);  // Row indices
    mwIndex *colPointers = mxGetJc(sparseArray); // Column pointers

    size_t numRows = static_cast<size_t>(mxGetM(sparseArray)); // Number of rows
    size_t numCols = static_cast<size_t>(mxGetN(sparseArray)); // Number of columns
    size_t nnz = static_cast<size_t>(colPointers[numCols]);    // Number of non-zero elements

    /* std::cout << "Dimenstion of the matrix: "
              << numRows << "x" << numCols
              << ", nnz = " << nnz << std::endl; */

    // Populate the sparse matrix
    std::vector<double> values(val, val + nnz);
    std::vector<ptrdiff_t> rows(rowIndices, rowIndices + nnz);
    std::vector<ptrdiff_t> cols(colPointers, colPointers + numCols + 1);

    sparse_matrix.mColPointers = std::move(cols);
    sparse_matrix.mRowIndices = std::move(rows);
    sparse_matrix.mValues = std::move(values);
    sparse_matrix.mNumCols = numCols;
    sparse_matrix.mNumRows = numRows;
    sparse_matrix.mNNZ = nnz;

    mxDestroyArray(sparseArray);
    matClose(matfile);
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
    std::string target_filename = config.filename + "ros2_A_0.mat";
    if (!read_mat(target_filename.c_str(), target_matrix)) {
        std::cerr << "Error reading the sparse matrix!" << std::endl;
        return 1;
    }

    // Read the reqired sparse matrices
    double simple_sparstiy_pattern_average_time = 0.0;
    double global_sparstiy_pattern_average_time = 0.0;
    double column_sparstiy_pattern_average_time = 0.0;
    double lflil_sparstiy_pattern_average_time = 0.0;
    for (int i = 1; i <= 160; i += 15)
    {
        std::string filename = config.filename + "ros2_A_" + std::to_string(i) + ".mat";
        csc_matrix<> test_matix;
        if (!read_mat(filename.c_str(), test_matix))
        {
            std::cerr << "Error reading the sparse matrix!" << std::endl;
            return 1;
        }

        // Run the benchmark for simple sparsity pattern
        double simple_sparstiy_pattern_time = 0.0;
        double simple_sparsity_pattern_sam_time = 0.0;
        int simple_sparsity_pattern_nnz = 0;
        for (int iter = 0; iter < config.iters; iter++)
        {
            sparsity_pattern<> sparsity_pattern;
            std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
            simple_sparsity_pattern(test_matix, sparsity_pattern);
            std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> time_per_iter = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
            simple_sparstiy_pattern_time += time_per_iter.count();
            simple_sparsity_pattern_nnz += sparsity_pattern.mNNZ;

            // Compute the SAM
            std::chrono::high_resolution_clock::time_point sam_start = std::chrono::high_resolution_clock::now();
            csc_matrix<> map;
            SAM(test_matix, target_matrix, sparsity_pattern, map);
            std::chrono::high_resolution_clock::time_point sam_end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> sam_time_per_iter = std::chrono::duration_cast<std::chrono::duration<double>>(sam_end - sam_start);
            simple_sparsity_pattern_sam_time += sam_time_per_iter.count();
        }
        simple_sparstiy_pattern_average_time += (simple_sparstiy_pattern_time + simple_sparsity_pattern_sam_time) / config.iters;


        // Run the benchmark for global threshold sparsity pattern
        double global_sparstiy_pattern_time = 0.0;
        double global_sparsity_pattern_sam_time = 0.0;
        int global_sparsity_pattern_nnz = 0;
        for (int iter = 0; iter < config.iters; iter++)
        {
            sparsity_pattern<> sparsity_pattern;
            std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
            sparsity_pattern_global_thresh(test_matix, 0.001, sparsity_pattern);
            std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> time_per_iter = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
            global_sparstiy_pattern_time += time_per_iter.count();
            global_sparsity_pattern_nnz += sparsity_pattern.mNNZ;

            // Compute the SAM
            std::chrono::high_resolution_clock::time_point sam_start = std::chrono::high_resolution_clock::now();
            csc_matrix<> map;
            SAM(test_matix, target_matrix, sparsity_pattern, map);
            std::chrono::high_resolution_clock::time_point sam_end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> sam_time_per_iter = std::chrono::duration_cast<std::chrono::duration<double>>(sam_end - sam_start);
            global_sparsity_pattern_sam_time += sam_time_per_iter.count();
        }
        global_sparstiy_pattern_average_time += (global_sparstiy_pattern_time + global_sparsity_pattern_sam_time) / config.iters;

        // Run the benchmark for column threshold sparsity pattern
        double column_sparstiy_pattern_time = 0.0;
        double column_sparsity_pattern_sam_time = 0.0;
        int column_sparsity_pattern_nnz = 0;
        for (int iter = 0; iter < config.iters; iter++)
        {
            sparsity_pattern<> sparsity_pattern;
            std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
            sparsity_pattern_col_thresh(test_matix, 0.7, sparsity_pattern);
            std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> time_per_iter = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
            column_sparstiy_pattern_time += time_per_iter.count();
            column_sparsity_pattern_nnz += sparsity_pattern.mNNZ;

            // Compute the SAM
            std::chrono::high_resolution_clock::time_point sam_start = std::chrono::high_resolution_clock::now();
            csc_matrix<> map;
            SAM(test_matix, target_matrix, sparsity_pattern, map);
            std::chrono::high_resolution_clock::time_point sam_end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> sam_time_per_iter = std::chrono::duration_cast<std::chrono::duration<double>>(sam_end - sam_start);
            column_sparsity_pattern_sam_time += sam_time_per_iter.count();
        }
        column_sparstiy_pattern_average_time += (column_sparstiy_pattern_time + column_sparsity_pattern_sam_time) / config.iters;

        // Run the benchmark for lfil threshold sparsity pattern
        double lfil_sparstiy_pattern_time = 0.0;
        double lfil_sparsity_pattern_sam_time = 0.0;
        int lfil_sparsity_pattern_nnz = 0;
        for (int iter = 0; iter < config.iters; iter++)
        {
            sparsity_pattern<> sparsity_pattern;
            std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
            sparsity_pattern_lfil_thresh(test_matix, 5, sparsity_pattern);
            std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> time_per_iter = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
            lfil_sparstiy_pattern_time += time_per_iter.count();
            lfil_sparsity_pattern_nnz += sparsity_pattern.mNNZ;

            // Compute the SAM
            std::chrono::high_resolution_clock::time_point sam_start = std::chrono::high_resolution_clock::now();
            csc_matrix<> map;
            SAM(test_matix, target_matrix, sparsity_pattern, map);
            std::chrono::high_resolution_clock::time_point sam_end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> sam_time_per_iter = std::chrono::duration_cast<std::chrono::duration<double>>(sam_end - sam_start);
            lfil_sparsity_pattern_sam_time += sam_time_per_iter.count();
        }
        lflil_sparstiy_pattern_average_time += (lfil_sparstiy_pattern_time + lfil_sparsity_pattern_sam_time) / config.iters;

        std::cout << "Simple sparsity pattern time, Matrix A_" << i << ": " << simple_sparstiy_pattern_time / config.iters << " seconds" << std::endl;
        std::cout << "Simple sparsity pattern SAM computation time, Matrix A_" << i << ": " << simple_sparsity_pattern_sam_time / config.iters << " seconds" << std::endl;
        std::cout << "Global threshold sparsity pattern time, Matrix A_" << i << ": " << global_sparstiy_pattern_time / config.iters << " seconds" << std::endl;
        std::cout << "Global threshold sparsity pattern SAM computation time, Matrix A_" << i << ": " << global_sparsity_pattern_sam_time / config.iters << " seconds" << std::endl;
        std::cout << "Column threshold sparsity pattern time, Matrix A_" << i << ": " << column_sparstiy_pattern_time / config.iters << " seconds" << std::endl;
        std::cout << "Column threshold sparsity pattern SAM computation time, Matrix A_" << i << ": " << column_sparsity_pattern_sam_time / config.iters << " seconds" << std::endl;
        std::cout << "lfil sparsity pattern time, Matrix A_" << i << ": " << lfil_sparstiy_pattern_time / config.iters << " seconds" << std::endl;
        std::cout << "lfil sparsity pattern SAM computation time, Matrix A_" << i << ": " << lfil_sparsity_pattern_sam_time / config.iters << " seconds" << std::endl;
        std::cout << std::endl;

        std::cout << "NNZ of matrix A_" << i << ": " << test_matix.mNNZ << std::endl;
        std::cout << "Simple sparsity pattern NNZ, Matrix A_" << i << ": " << simple_sparsity_pattern_nnz / config.iters << std::endl;
        std::cout << "Global threshold sparsity pattern NNZ, Matrix A_" << i << ": " << global_sparsity_pattern_nnz / config.iters << std::endl;
        std::cout << "Column threshold sparsity pattern NNZ, Matrix A_" << i << ": " << column_sparsity_pattern_nnz / config.iters << std::endl;
        std::cout << "lfil sparsity pattern NNZ, Matrix A_" << i << ": " << lfil_sparsity_pattern_nnz / config.iters << std::endl;
        std::cout << std::endl;
    }

    std::cout << "Simple sparsity pattern (average computation time): " << simple_sparstiy_pattern_average_time / 10 << " seconds" << std::endl;
    std::cout << "Global threshold sparsity pattern (average computation time): " << global_sparstiy_pattern_average_time / 10 << " seconds" << std::endl;
    std::cout << "Column threshold sparsity pattern (average computation time): " << column_sparstiy_pattern_average_time / 10 << " seconds" << std::endl;
    std::cout << "lfil sparsity pattern (average computation time): " << lflil_sparstiy_pattern_average_time / 10 << " seconds" << std::endl;
    std::cout << std::endl;

    return 0;
}