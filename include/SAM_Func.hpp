#pragma once

#include <vector>
#include <mat.h>
#include <iostream>
#include <set>

class csc_matrix
{
public:
    csc_matrix() = default;
    csc_matrix(const std::vector<double> &values, const std::vector<size_t> &rowIndices, const std::vector<size_t> &colPointers, std::vector<size_t> &NNZPerCol, size_t numCols, size_t numRows, size_t nnz)
        : mValues(values), mRowIndices(rowIndices), mColPointers(colPointers), mNNZPerCol(NNZPerCol), mNumCols(numCols), mNumRows(numRows), mNNZ(nnz) {}

    csc_matrix(const csc_matrix &other)
    {
        mValues = other.getValues();
        mRowIndices = other.getRowIndices();
        mColPointers = other.getColPointers();
        mNNZPerCol = other.getNNZPerCol();
        mNumCols = other.mNumCols;
        mNumRows = other.mNumRows;
        mNNZ = other.mNNZ;
    }

    // Getters
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
    std::vector<double> mValues;
    std::vector<size_t> mRowIndices;
    std::vector<size_t> mColPointers;
    std::vector<size_t> mNNZPerCol;
    size_t mNumCols;
    size_t mNumRows;
    size_t mNNZ;
};

void simple_sparsity_pattern(const csc_matrix& A, csc_matrix& S)
{
    S.setRowIndices(A.getRowIndices());
    S.setColPointers(A.getColPointers());
    S.setNumCols(A.getNumCols());
    S.setNumRows(A.getNumRows());
    S.setNNZ(A.getNNZ());
    S.setNNZPerCol(A.getNNZPerCol());

    std::vector<double> values(A.getNNZ(), 1);
    S.setValues(std::cref(values));
}

void extractSubmatrix(const std::vector<size_t>& submatrixColIndices, std::vector<size_t>& submatrixRowIndices, std::vector<size_t>& submatrixRowPointers, const std::vector<size_t>& submatrixColPointers, const csc_matrix& S, int& maxSk, int& maxRk) {
    // Preprocessing the sparsity pattern
    const size_t numCols = S.getNumCols();
    const std::vector<size_t> nnzPerCol = S.getNNZPerCol();

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

        // Iterate over all the non zero indices in the current columna and add the nonzero indices of the respecticve columns to the set
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


csc_matrix SAM(const csc_matrix& source, const csc_matrix& target, const csc_matrix& S)
{
    csc_matrix MM;

    // Construct the submatrix information for each column
    std::vector<size_t> submatrixColIndices(S.getRowIndices());
    std::vector<size_t> submatrixColPointers(S.getColPointers());
    std::vector<size_t> submatrixRowIndices;
    std::vector<size_t> submatrixRowPointers{0};

    // Keep track of the maximum number of non zeros in the row and columns of the submatrix
    int maxSk = 0;
    int maxRk = 0;
    extractSubmatrix(submatrixColIndices, submatrixRowIndices, submatrixRowPointers, submatrixColPointers, S, maxSk, maxRk);

    return MM;
}