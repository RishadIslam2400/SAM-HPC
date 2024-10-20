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
void sparsity_pattern_global_thresh(const csc_matrix &A, double thresh, csc_matrix &S) {
    const std::vector<double> inputVals = A.getValuesRef();
    const std::vector<size_t> inputRowIndices = A.getRowIndicesRef();
    const std::vector<size_t> inputColPointers = A.getColPointersRef();
    const size_t numCols = A.getNumCols();
    const size_t numRows = A.getNumRows();
    const size_t nnz = A.getNNZ();

    // TODO: Parallelize
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
}