#include "cscMatrix.hpp"
#include "config.hpp"
#include "sparsityPattern.hpp"

#include <mat.h>
#include <matrix.h>
#include <iostream>
#include <vector>
#include <functional>
#include <chrono>

// This benchmark will compare the performance between all the different sparsity pattern choices
// We will compare the different sparsity pattern choices as a whole
// and then look into each spasity pattern to see which part of the code is taking the most time

bool read_mat(const char *filename, csc_matrix<double, ptrdiff_t> &sparse_matrix) {
    // Read the .mat file
    MATFile *matfile = matOpen(filename, "r");
    if (!matfile) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return false;
    }

    // Get all variable names in the file
    int num_variables;
    const char **dir = (const char **)matGetDir(matfile, &num_variables);
    if (!dir) {
        std::cerr << "Error getting variables" << std::endl;
        matClose(matfile);
        return false;
    }

    // Read the sparse matrix variable
    mxArray *sparseArray = nullptr;

    // Loop through the variables to find a sparse matrix
    for (int i = 0; i < num_variables; ++i) {
        mxArray *var = matGetVariable(matfile, dir[i]);
        if (!var) {
            std::cerr << "Error getting variable: " << dir[i] << std::endl;
            mxFree(dir);
            matClose(matfile);
            return false;
        }

        // Check if the variable is a sparse matrix
        if (mxIsSparse(var)) {
            sparseArray = var; // Found the sparse matrix
            break;             // Stop searching
        }

        mxDestroyArray(var); // Free memory for non-sparse variables
    }
    mxFree(dir);

    // Ensure the variable is a sparse matrix
    if (!sparseArray) {
        std::cerr << "No sparse matrix found in the file!" << std::endl;
        matClose(matfile);
        return false;
    }

    // Get the data - Matlab stores sparse matrix in csc format
    double *val = mxGetPr(sparseArray);                        // Non-zero elements
    mwIndex *rowIndices = mxGetIr(sparseArray);                // Row indices
    mwIndex *colPointers = mxGetJc(sparseArray);               // Column pointers

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

int main(int argc, char **argv) {
    // get the configuration, print it
    config_t config;
    parseargs(argc, argv, config);
    config.print();

    // Read the reqired sparse matrices
    csc_matrix<double, ptrdiff_t> test_matix;
    if (!read_mat(config.filename.c_str(), test_matix)) {
        std::cerr << "Error reading the sparse matrix!" << std::endl;
        return 1;
    }

    // Run the benchmark for global sparsity pattern
    double function_times[3] = {0.0, 0.0, 0.0};
    for (int iter = 0; iter < config.iters; iter++) {
        sparsity_pattern<> test_sparsity_pattern;
        sparsity_pattern_global_thresh(test_matix, 0.001, test_sparsity_pattern, function_times);
    }

    std::cout << "Diagonal scaling time: " << function_times[0] / config.iters << " seconds" << std::endl;
    std::cout << "Sparsity pattern construction time: " << function_times[1] / config.iters << " seconds" << std::endl;
    std::cout << "Total time: " << function_times[2] / config.iters << " seconds" << std::endl;

    
}