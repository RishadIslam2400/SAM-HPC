#include <iostream>
#include <vector>
#include <mat.h>
#include <matrix.h>
#include <string>

#include "sparsityPattern.hpp"
#include "samFunc.hpp"
#include "cscMatrix.hpp"

// Function to read .mat file directly using MATLAB API
bool read_mat(const char *filename, csc_matrix<double, ptrdiff_t>& sparse_matrix) {
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

    std::cout << "Dimenstion of the matrix: "
              << numRows << "x" << numCols
              << ", nnz = " << nnz << std::endl;


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

int main() {
    csc_matrix<> A;
    sparsity_pattern<> S;

    const std::string fileName = "/home/rishad/SAM-HPC/data/matrix_1.mat";
    if (!read_mat(fileName.c_str(), A)) {
        std::cerr << "Error reading matrix file" << std::endl;
        return -1;
    }

    // sparsity_pattern_global_thresh(A, 0.001, S);
    // sparsity_pattern_col_thresh(A, 0.9, S);
    sparsity_pattern_lfil_thresh(A, 3, S);

    S.printMatrix();
    // S.printColumn(1);
}