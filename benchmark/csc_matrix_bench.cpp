#include "cscMatrix.hpp"
#include "config.hpp"

#include <mat.h>
#include <matrix.h>
#include <iostream>
#include <vector>
#include <functional>
#include <chrono>

// This benchmark will test the performance of differnt functions within the cscMatrix class.
// The perforance will be compared between matrices of differnt sizes.
// All the implementations within the csc matrix class are sequential.

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

    // Benchmarking the read_mat function - time taken to read each matrix
    /* std::chrono::high_resolution_clock::time_point read_start = std::chrono::high_resolution_clock::now();
    for (int iter = 0; iter < config.iters; iter++) {
        csc_matrix<double, ptrdiff_t> test_matrix;
        if (!read_mat(config.filename.c_str(), test_matrix)) {
            std::cerr << "Error reading matrix file" << std::endl;
            return -1;
        }
    }
    std::chrono::high_resolution_clock::time_point read_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> read_time = std::chrono::duration_cast<std::chrono::duration<double>>(read_end - read_start);
    std::cout << "Time taken to read matrices: " << read_time.count() / (config.iters) << " seconds" << std::endl;

    std::cout << std::endl; */

    // Populate the test matrix
    csc_matrix<double, ptrdiff_t> test_matrix;
    if (!read_mat(config.filename.c_str(), test_matrix)) {
        std::cerr << "Error reading matrix file" << std::endl;
        return -1;
    }

    // Benchmarking the extract diagonal and diagonal scaling functions
    // double time_spans[4] = {0.0, 0.0, 0.0, 0.0};
    // std::chrono::high_resolution_clock::time_point scaling_start = std::chrono::high_resolution_clock::now();
    for (int iter = 0; iter < config.iters; iter++) {
        // std::vector<std::chrono::high_resolution_clock::time_point> times(6);
        // times[0] = std::chrono::high_resolution_clock::now();
        // Perform diagonal scaling of the matrix values
        auto scaledValues = test_matrix.diagonalScale();
        // times[5] = std::chrono::high_resolution_clock::now();

        // time_spans[0] += std::chrono::duration_cast<std::chrono::duration<double>>(times[5] - times[0]).count();
        // time_spans[1] += std::chrono::duration_cast<std::chrono::duration<double>>(times[1] - times[0]).count();
        // time_spans[2] += std::chrono::duration_cast<std::chrono::duration<double>>(times[3] - times[2]).count();
        // time_spans[3] += std::chrono::duration_cast<std::chrono::duration<double>>(times[4] - times[3]).count();
    }

    // std::chrono::high_resolution_clock::time_point scaling_end = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> scaling_time_span = std::chrono::duration_cast<std::chrono::duration<double>>(scaling_end - scaling_start);
    // std::cout << "Time taken to perfrom complete diagonal scaling: " << time_spans[0] / (config.iters) << " seconds" << std::endl;
    // std::cout << "Time taken to extract the diagonal: " << time_spans[1] / (config.iters) << " seconds" << std::endl;
    // std::cout << "Time taken to post multiply the matrix: " << time_spans[2] / (config.iters) << " seconds" << std::endl;
    // std::cout << "Time taken to pre multiply the matrix: " << time_spans[3] / (config.iters) << " seconds" << std::endl;

    std::cout << std::endl;

    return 0;
}