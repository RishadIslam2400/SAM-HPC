#pragma once

#include <vector>
#include <mat.h>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <iomanip>

// Representation of the sparse matrices
// It uses compressed sparse column (CSC) storage format
class csc_matrix
{
public:
    csc_matrix() = default;
    
    // Constructor
    csc_matrix(const std::vector<double> &values, const std::vector<size_t> &rowIndices, const std::vector<size_t> &colPointers, std::vector<size_t> &NNZPerCol, size_t numCols, size_t numRows, size_t nnz)
        : mValues(values), mRowIndices(rowIndices), mColPointers(colPointers), mNNZPerCol(NNZPerCol), mNumCols(numCols), mNumRows(numRows), mNNZ(nnz) {}

    // Copy constructor
    csc_matrix(const csc_matrix &other)
        : mValues(other.mValues), mRowIndices(other.mRowIndices), mColPointers(other.mColPointers),
          mNNZPerCol(other.mNNZPerCol), mNumCols(other.mNumCols), mNumRows(other.mNumRows), mNNZ(other.mNNZ) {}

    // Getters - the vectors are const references 
    const std::vector<double> &getValues() const { return mValues; }
    const std::vector<size_t> &getRowIndices() const { return mRowIndices; }
    const std::vector<size_t> &getColPointers() const { return mColPointers; }
    const std::vector<size_t> &getNNZPerCol() const { return mNNZPerCol; }
    size_t getNumCols() const { return mNumCols; }
    size_t getNumRows() const { return mNumRows; }
    size_t getNNZ() const { return mNNZ; }

    // Setters
    void setValues(const std::vector<double> &values) { mValues = values; }
    void setRowIndices(const std::vector<size_t> &rowIndices) { mRowIndices = rowIndices; }
    void setColPointers(const std::vector<size_t> &colPointers) { mColPointers = colPointers; }
    void setNNZPerCol(const std::vector<size_t> &NNZPerCol) { mNNZPerCol = NNZPerCol; }
    void setNumCols(size_t numCols) { mNumCols = numCols; }
    void setNumRows(size_t numRows) { mNumRows = numRows; }
    void setNNZ(size_t nnz) { mNNZ = nnz; }

    void printMatrix() const
    {
        for (size_t col = 0; col < mNumCols; ++col)
        {
            for (size_t idx = mColPointers[col]; idx < mColPointers[col + 1]; ++idx)
            {
                std::cout << "(" << mRowIndices[idx] + 1 << ", " << col + 1 << ") = " << mValues[idx] << std::endl;
            }
        }
    }

private:
    std::vector<double> mValues;      // Non-zero elements
    std::vector<size_t> mRowIndices;  // Row indices
    std::vector<size_t> mColPointers; // Pointer to the start of each column in the row indices array
    std::vector<size_t> mNNZPerCol;   // Number of non zero elements in each column
    size_t mNumCols;                  // total number of columns
    size_t mNumRows;                  // total number of rows
    size_t mNNZ;                      // total number of non-zero elements
};

// Function to generate a simple sparsity pattern
// The sparsity pattern is same as the input matrix
void simple_sparsity_pattern(const csc_matrix& A, csc_matrix& S)
{
    S.setRowIndices(A.getRowIndices());
    S.setColPointers(A.getColPointers());
    S.setNumCols(A.getNumCols());
    S.setNumRows(A.getNumRows());
    S.setNNZ(A.getNNZ());
    S.setNNZPerCol(A.getNNZPerCol());

    // Instead of the values of the original matrices, the non zero elements are set to 1
    std::vector<double> values(A.getNNZ(), 1);
    S.setValues(std::move(values));
}

// Function to extract the submatrix information from the sparsity pattern and the source matrix
void extractSubmatrixInfo(const std::vector<size_t>& submatrixColIndices, std::vector<size_t>& submatrixRowIndices, std::vector<size_t>& submatrixRowPointers, const std::vector<size_t>& submatrixColPointers, const csc_matrix& S, int& maxSk, int& maxRk) {
    // Preprocessing the sparsity pattern
    const size_t numCols = S.getNumCols();
    const std::vector<size_t>& nnzPerCol = S.getNNZPerCol();

    // Reserve memory
    submatrixRowPointers.reserve(numCols + 1);

    // Naive implementation
    // For each column find the nonzero indices
    // For each column corresponding to the nonzero indices get a set union of non zero indices for the row indices 
    for (size_t j = 0; j < numCols; j++) {
        // Static allocation is difficult since we do not know the size beforehand
        // Use set to keep track of the unique row indices
        std::set<size_t> rowIndicesPerCol;
        const auto start = submatrixColIndices.begin() + submatrixColPointers[j];
        const auto end = submatrixColIndices.begin() + submatrixColPointers[j + 1];

        // Update the maxSk value for the current column
        maxSk = std::max(maxSk, static_cast<int>(nnzPerCol[j]));

        // Iterate over all the non zero indices in the current column and add the nonzero indices of the respecticve columns to the set
        for (auto it = start; it != end; ++it) {
            const size_t col = *it;
            auto colStart = submatrixColIndices.begin() + submatrixColPointers[col];
            auto colEnd = submatrixColIndices.begin() + submatrixColPointers[col + 1];
            rowIndicesPerCol.insert(colStart, colEnd);
        }

        // Update the maxRk value for the current column
        maxRk = std::max(maxRk, static_cast<int>(rowIndicesPerCol.size()));

        // Insert into the row indices for the submatrix for the current column
        submatrixRowIndices.insert(submatrixRowIndices.end(), rowIndicesPerCol.begin(), rowIndicesPerCol.end());

        // Keep track of the number of non zero elements in the submatrix
        submatrixRowPointers.push_back(submatrixRowIndices.size());
    }
}


// The SAM algorithm
csc_matrix SAM(const csc_matrix& source, const csc_matrix& target, const csc_matrix& S)
{
    csc_matrix MM;

    // Construct the submatrix information for each column
    const std::vector<size_t> submatrixColIndices(S.getRowIndices());
    const std::vector<size_t> submatrixColPointers(S.getColPointers());
    std::vector<size_t> submatrixRowIndices;
    std::vector<size_t> submatrixRowPointers{0};

    // Get other information from the sparsity pattern
    const size_t numCols = S.getNumCols();
    const std::vector<size_t> nnzPerCol = S.getNNZPerCol();

    // Keep track of the maximum number of non zeros in the row and columns of the submatrix
    int maxSk = 0;
    int maxRk = 0;
    extractSubmatrixInfo(submatrixColIndices, submatrixRowIndices, submatrixRowPointers, submatrixColPointers, S, maxSk, maxRk);

    // Get the values and row indices from the source matrix
    const std::vector<double>& sourceMatrixValues(source.getValues());
    const std::vector<size_t>& sourceMatrixRowIndices(source.getRowIndices());


    // Extract the submatrix from source matrix using the submatrix information
    for (size_t i = 0; i < source.getNumCols(); i++) {
        // Find the column and row indices of the submatrix for the current column
        auto colStart = submatrixColIndices.begin() + submatrixColPointers[i];
        size_t colDim = submatrixColPointers[i + 1] - submatrixColPointers[i];

        auto rowStart = submatrixRowIndices.begin() + submatrixRowPointers[i];
        size_t rowDim = submatrixRowPointers[i + 1] - submatrixRowPointers[i];

        // Initialize the submatrix with fixed size
        // The submatrices are in row major order
        // std::vector<std::vector<double>> submatrix(rowDim, std::vector<double>(colDim));

        // This submatrix will be in column major order
        std::vector<std::vector<double>> submatrix(colDim, std::vector<double>(rowDim));

        // Map row indices to their position in submatrix
        std::unordered_map<size_t, size_t> rowIndexMap;
        for (size_t row = 0; row < rowDim; ++row, ++rowStart) {
            rowIndexMap[*rowStart] = row;
        }

        // Iterate over the column indices to build the submatrix
        for (auto j = 0; j < colDim; ++j, ++colStart) {
            // // Get the iterators for the non-zero values of the column
            const size_t col = *colStart;
            auto start = sourceMatrixValues.begin() + submatrixColPointers[col];
            auto end = sourceMatrixValues.begin() + submatrixColPointers[col + 1];
            const size_t nnz = submatrixColPointers[col + 1] - submatrixColPointers[col];

            // Get the row indices for the column
            auto rowIndexStart = sourceMatrixRowIndices.begin() + submatrixColPointers[col];

            // Check the row index for the current column and insert the non zero values in the exact position in the submatrix
            for (auto k = start; k != end; ++k, ++rowIndexStart) {
                // accessing row major submatrix
                // submatrix[rowIndexMap[*rowIndexStart]][j] = *k;

                // accessing column major submatrix
                submatrix[j][rowIndexMap[*rowIndexStart]] = *k;
            }
        }

        // Solve the submatrix using QR factorization with Householder Transformations

        // Extract the required required information from the target matrix for solving the LS problem
        const std::vector<double>& targetMatrixValues(target.getValues());
        const std::vector<size_t>& targetMatrixRowIndices(target.getRowIndices());
        const std::vector<size_t>& targetMatrixColPointers(target.getColPointers());

        // Extract current column from the target matrix
        auto targetColStart = targetMatrixValues.begin() + targetMatrixColPointers[i];
        auto targetColEnd = targetMatrixValues.begin() + targetMatrixColPointers[i + 1];

        // Initialize the RHS vector maintaining the dimension of the submatrix
        std::vector<double> rhs(rowDim);

        // Generate the RHS vector
        auto targetRowIndexStart = targetMatrixRowIndices.begin() + targetMatrixColPointers[i];
        for (auto k = targetColStart; k != targetColEnd; ++k, ++targetRowIndexStart) {
            rhs[rowIndexMap[*targetRowIndexStart]] = *k;
        }

        // Initialize the map vector to store the solution
        std::vector<double> mapColumn(colDim);

        // Solve the submatrix using QR factorization with Householder Transformations
        // Pass the arguments by value to avoid modifying the original matrices
        // qrHouseholder(submatrix, rhs, mapColumn, colDim, rowDim);

        // Debug Print: Print the submatrix for the current column
        std::cout << "Submatrix for column " << i << ":\n";
        for (int i = 0; i < rowDim; ++i) {
            for (int j = 0; j < colDim; ++j) {
                std::cout << std::setw(10) << std::setprecision(10) << submatrix[j][i] << " ";
            }
            std::cout << "\n";
        }
        std::cout << "\n";
    }

    return MM;
}