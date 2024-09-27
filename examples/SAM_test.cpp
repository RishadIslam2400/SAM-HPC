#include <mat.h>
#include <matrix.h>
#include <iostream>
#include <vector>
#include <tuple>
#include <functional>

#include "SAM_Func.hpp"

// Function to read .mat file directly using MATLAB API
bool read_mat(const char *filename, csc_matrix& sparse_matrix)
{
    // Read the .mat file
    MATFile *matfile = matOpen(filename, "r");
    if (!matfile)
    {
        std::cerr << "Error opening file" << std::endl;
        return false;
    }

    // Read the sparse matrix variable
    mxArray  *a = matGetVariable(matfile, "Jac");
    if (!a) {
        std::cerr << "Error getting variable" << std::endl;
        matClose(matfile);
        return false;
    }

    // Ensure the variable is a sparse matrix
    if (!mxIsSparse(a)) {
        std::cerr << "The variable is not a sparse matrix!" << std::endl;
        mxDestroyArray(a);
        matClose(matfile);
        return false;
    }

    // Get the data
    double *val = mxGetPr(a);                               // Non-zero elements
    mwIndex *rowIndices = mxGetIr(a);                       // Row indices
    mwIndex *colPointers = mxGetJc(a);                      // Column pointers
    size_t numRows = static_cast<size_t>(mxGetM(a));        // Number of rows
    size_t numCols = static_cast<size_t>(mxGetN(a));        // Number of columns
    size_t nnz = static_cast<size_t>(colPointers[numCols]); // Number of non-zero elements

    // Calculate the number of non-zero elements per column
    std::vector<size_t> nnzPerCol(numCols);
    for (size_t col = 0; col < numCols; ++col) {
        nnzPerCol[col] = colPointers[col + 1] - colPointers[col];
    }

    // Populate the sparse matrix
    std::vector<double> values(val, val + nnz);
    std::vector<size_t> rows(rowIndices, rowIndices + nnz);
    std::vector<size_t> cols(colPointers, colPointers + numCols + 1);
    
    sparse_matrix.setValues(values);
    sparse_matrix.setRowIndices(rows);
    sparse_matrix.setColPointers(cols);
    sparse_matrix.setNNZPerCol(nnzPerCol);
    sparse_matrix.setNumCols(numCols);
    sparse_matrix.setNumRows(numRows);
    sparse_matrix.setNNZ(nnz);

    mxDestroyArray(a);
    matClose(matfile);
    return true;
}

int main()
{
    // Seqeunce of matrices
    std::vector<csc_matrix> sequence;

    // Number of matrices. Currently specified as per the requirement of the application
    unsigned int numMatrices = 70;

    // Read the sequence of the matrices from .mat files and store them in the vector
    for (int i = 1 ; i <= numMatrices; i++) {
        std::string fileName = "/home/rishad/SAM-HPC/data/matrix_" + std::to_string(i) + ".mat";
        csc_matrix temp;
        if (read_mat(fileName.c_str(), temp)) {
            sequence.emplace_back(std::move(temp));
        }
        else {
            std::cerr << "Error reading matrix file" << std::endl;
            return -1;
        }
    }

    // The initial matrix is the target matrix
    csc_matrix A0(std::cref(sequence[0]));    

    // SAM function to map other matrices in the sequence to the target matrix
    for (int i = 1; i < 2; i++) {
        csc_matrix Ak(std::cref(sequence[i])); // Current source matrix
        csc_matrix Sk;                         // Sparsity pattern for the SAM
        csc_matrix Mk;                         // Actual values of the map
        
        // Compute the sparsity pattern for SAM for the current source matrix
        simple_sparsity_pattern(A0, Sk);
        
        // Compute the map
        Mk = SAM(A0, Ak, Sk);

        // Print the map matrix
        std::cout << "Mk: " << std::endl;
        Mk.printMatrix();
        std::cout << std::endl;
    }

    return 0;    
}