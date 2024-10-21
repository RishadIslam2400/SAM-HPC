#include "samFunc.hpp"


// Function to generate a simple sparsity pattern
// The sparsity pattern is same as the input matrix
void simple_sparsity_pattern(const csc_matrix &A, csc_matrix &S)
{
    S.setRowIndices(A.getRowIndicesRef());
    S.setColPointers(A.getColPointersRef());
    S.setNumCols(A.getNumCols());
    S.setNumRows(A.getNumRows());
    S.setNNZ(A.getNNZ());

    // Instead of the values of the original matrices, the non zero elements are set to 1
    std::vector<double> values(A.getNNZ(), 1);
    S.setValues(std::move(values));
}

// Generates a sparsity pattern based global threshold parameter
// Sparsifies the input matrix A based on the global threshold
// Matrix entries less than the global threshold are set to 0
// For better approximation, level 2 neighbars of the sparse matrix is taken (power of 2)
// Post sparsification, for further optimization
void sparsity_pattern_global_thresh(const csc_matrix &A,const double thresh, csc_matrix &S) {
    std::vector<double> inputVals{A.getValuesCopy()};
    std::vector<size_t> inputRowIndices{A.getRowIndicesCopy()};
    std::vector<size_t> inputColPointers{A.getColPointersCopy()};
    const size_t numCols = A.getNumCols();
    const size_t numRows = A.getNumRows();
    const size_t nnz = A.getNNZ();

    // TODO: Parallelize
    // TODO: Integrate these matrix operations with the class definition
    // Extract the diagonal entries from A for diagonal scaling
    // Take the absolute value of the diagonal and compute the inverse of the square root
    std::vector<double> diag(numCols, 0.0);
    for (size_t j = 0; j < numCols; ++j) {
        // Iterate over the non-zeros of the current column
        size_t colStart = inputColPointers[j];
        size_t colEnd = inputColPointers[j + 1];
        for (size_t i = colStart; i < colEnd; ++i) {
            if (inputRowIndices[i] == j) {
                double val = std::abs(inputVals[i]);
                val = 1.0 / std::sqrt(val);
                diag[j] = val;
                break;
            }
        }
        // If the diagonal entry is 0 then replace it with 1
        if (diag[j] == static_cast<double>(0)) {
            diag[j] = static_cast<double>(1);
        }
    }

    // TODO: scale the matrix (value array) with the diagonal array, A_norm = D^-1/2 * A * D^-1/2
    // TODO: Parallelize
    // Post multiplying a matrix by a diagonal matrix,
    // [a1] = d1 * [a1] (multiplying each diagonal element with the corresponding column in the matrix)
    for (size_t j = 0; j < numCols; ++j) {
        size_t colStart = inputColPointers[j];
        size_t colEnd = inputColPointers[j + 1];
        for (size_t i = colStart; i < colEnd; ++i) {
            inputVals[i] *= diag[j];
        }
    }

    // Pre multiplying a matrix by a diagonal matrix,
    // [a1]^T = d1 * [a1]^T (multiplying each diagonal element with the corresponding row in the matrix)
    for (size_t j = 0; j < numCols; ++j) {
        size_t colStart = inputColPointers[j];
        size_t colEnd = inputColPointers[j + 1];
        for (size_t i = colStart; i < colEnd; ++i) {
            size_t idx = inputRowIndices[i];
            inputVals[i] *= diag[idx]; // instead of iterating over the row, multiply with corresponding diagonal element
        }
    }

    // Filter the values with the threshold with the input values array
    // Genrate the val, rowIndices and colPointers for the sparsity pattern
    // TODO: Use arrays instead of vectors for better performance
    // TODO: change the type of the sparsity pattern values by templating
    std::vector<double> sparsityValues;
    std::vector<size_t> sparsityRowIndices;
    std::vector<size_t> sparsityColPointers;
    sparsityValues.reserve(nnz);
    sparsityRowIndices.reserve(nnz);
    sparsityColPointers.reserve(numCols + 1);
    sparsityColPointers.push_back(0);

    for (size_t j = 0; j < numCols; ++j) {
        const size_t colStart = inputColPointers[j];
        const size_t colEnd = inputColPointers[j + 1];
        unsigned int nnzCount{0};
        for (size_t i = colStart; i < colEnd; ++i) {
            if (std::abs(inputVals[i]) > thresh) {
                // TODO: Optmize the values array by initilizing later
                sparsityValues.push_back(1.0);
                sparsityRowIndices.push_back(inputRowIndices[i]);
                ++nnzCount;
            }
        }
        sparsityColPointers.push_back(nnzCount);
    }

    const size_t sparsityNNZ = sparsityRowIndices.size();
    const size_t sparsityNumCols = numCols;
    const size_t sparsityNumRows = numRows;

    // TODO: Resize the vectors if too small compared to the input vector

    // Construct the sparsity pattern matrix
    S.setNNZ(sparsityNNZ);
    S.setNumCols(sparsityNumCols);
    S.setNumRows(sparsityNumRows);
    S.setRowIndices(std::move(sparsityRowIndices));
    S.setColPointers(std::move(sparsityColPointers));
    S.setValues(std::move(sparsityValues));

    // TODO: Matrix matrix multiplication of S for level 2 neighbors
}

// Naive implementation
void sparsity_pattern_col_thresh(const csc_matrix &A, const double tau, csc_matrix &S) {
    std::vector<double> inputVals{A.getValuesCopy()};
    std::vector<size_t> inputRowIndices{A.getRowIndicesCopy()};
    std::vector<size_t> inputColPointers{A.getColPointersCopy()};
    const size_t numCols = A.getNumCols();
    const size_t numRows = A.getNumRows();
    const size_t nnz = A.getNNZ();

    // Information to construct the sparsity pattern
    std::vector<double> sparsityValues;
    std::vector<size_t> sparsityRowIndices;
    std::vector<size_t> sparsityColPointers;
    sparsityValues.reserve(nnz);
    sparsityRowIndices.reserve(nnz);
    sparsityColPointers.reserve(numCols + 1);
    sparsityColPointers.push_back(0);

    // Sparsify each column for the sparsity pattern S
    for (size_t j = 0; j < numCols; ++j) {
        const size_t colStart = inputColPointers[j];
        const size_t colEnd = inputColPointers[j + 1];
        unsigned int nnzCount{0};

        // Find the maximum absolute value in the current column
        double maxVal{0.0};
        for (size_t i = colStart; i < colEnd; ++i) {
            if (std::abs(inputVals[i]) > maxVal) {
                maxVal = std::abs(inputVals[i]);
            }
        }

        // Find the threshold for the current column
        const double colThresh{(1.0 - tau) * maxVal};

        // Sparsify the column and construct the sparsity pattern
        for (size_t i = colStart; i < colEnd; ++i) {
            // Keep the values greater than the threshold or the diagonal entry
            if (std::abs(inputVals[i]) > colThresh || inputRowIndices[i] == j) {
                sparsityValues.push_back(1.0);
                sparsityRowIndices.push_back(inputRowIndices[i]);
                ++nnzCount;
            }
        }
        sparsityColPointers.push_back(nnzCount);
    }

    const size_t sparsityNNZ = sparsityRowIndices.size();
    const size_t sparsityNumCols = numCols;
    const size_t sparsityNumRows = numRows;

    // Debug Print
    /* std::cout << "sparsityNNZ: " << sparsityNNZ << std::endl;
    std::cout << "rowIndices: " << std::endl;
    for (size_t i = 0; i < sparsityNNZ; ++i) {
        std::cout << sparsityRowIndices[i] << " ";
    }
    std::cout << std::endl;
    std::cout << "colPointers: " << std::endl;
    for (size_t i = 0; i < sparsityNumCols + 1; ++i) {
        std::cout << sparsityColPointers[i] << " ";
    }
    std::cout << std::endl; */

    // TODO: sparse matrix matrix multiplication of S for level 2 neighbors

    // Construct the sparsity pattern matrix
    S.setNNZ(sparsityNNZ);
    S.setNumCols(sparsityNumCols);
    S.setNumRows(sparsityNumRows);
    S.setRowIndices(std::move(sparsityRowIndices));
    S.setColPointers(std::move(sparsityColPointers));
    S.setValues(std::move(sparsityValues));
}