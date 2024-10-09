#include <mat.h>
#include <matrix.h>
#include <iostream>
#include <vector>
#include <tuple>
#include <functional>
#include <chrono>

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

    // Get all variable names in the file
    int num_variables;
    const char **dir = (const char **)matGetDir(matfile, &num_variables);
    if (!dir)
    {
        std::cerr << "Error getting variables directory" << std::endl;
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

    // Get the data
    double *val = mxGetPr(sparseArray);                     // Non-zero elements
    mwIndex *rowIndices = mxGetIr(sparseArray);             // Row indices
    mwIndex *colPointers = mxGetJc(sparseArray);                      // Column pointers
    size_t numRows = static_cast<size_t>(mxGetM(sparseArray));        // Number of rows
    size_t numCols = static_cast<size_t>(mxGetN(sparseArray));        // Number of columns
    size_t nnz = static_cast<size_t>(colPointers[numCols]); // Number of non-zero elements

    std::cout << "Number of rows: " << numRows << std::endl;
    std::cout << "Number of columns: " << numCols << std::endl;

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

    mxDestroyArray(sparseArray);
    matClose(matfile);
    return true;
}

int main()
{
    // Seqeunce of matrices
    std::vector<csc_matrix> sequence;

    // Number of matrices. Currently specified as per the requirement of the application
    unsigned int numMatrices = 3;

    // Read the sequence of the matrices from .mat files and store them in the vector
    for (int i = 1 ; i <= numMatrices; i++) {
        std::string fileName = "/home/rishad/SAM-HPC/data_13x21/ros2_A_" + std::to_string(i-1) + ".mat";
        csc_matrix temp;
        if (read_mat(fileName.c_str(), std::ref(temp))) {
            sequence.emplace_back(std::move(temp));
        }
        else {
            std::cerr << "Error reading matrix file" << std::endl;
            return -1;
        }
    }

    // The initial matrix is the target matrix
    csc_matrix A0(std::cref(sequence[0]));    

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

    // SAM function to map other matrices in the sequence to the target matrix
    for (int i = 1; i < numMatrices; i++) {
        csc_matrix Ak(std::cref(sequence[i])); // Current source matrix
        csc_matrix Sk;                         // Sparsity pattern for the SAM
        csc_matrix Mk;                         // Actual values of the map
        
        // Compute the sparsity pattern for SAM for the current source matrix
        simple_sparsity_pattern(A0, Sk);
        
        // Compute the map
        Mk = SAM(A0, Ak, Sk);

        // Print the map matrix
        /* std::cout << "Mk: " << std::endl;
        Mk.printMatrix();
        std::cout << std::endl; */
    }

    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    std::cout << "Time taken for each SAM computation: " << time_span.count() / (numMatrices - 1) << " seconds" << std::endl;

    return 0;    
}