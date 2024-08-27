#pragma once

#include <vector>
#include <mat.h>
#include <iostream>

// Class for the sparse matrix
class csc_matrix
{
public:
    csc_matrix() = default;
    csc_matrix(const std::vector<double> &values, const std::vector<size_t> &rowIndices, const std::vector<size_t> &colPointers, std::vector<size_t> &NNZPerCol, size_t numCols, size_t numRows, size_t nnz)
        : mValues(values), mRowIndices(rowIndices), mColPointers(colPointers), mNNZPerCol(NNZPerCol), mNumCols(numCols), mNumRows(numRows), mNNZ(nnz) {}

    // Copy constructor
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

// Sparsity pattern is same as the provided matrix
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


// Function to compute the maps. Contains both preprocessing and actual computation
csc_matrix SAM(const csc_matrix& source, const csc_matrix& target, const csc_matrix& S)
{
    // Preprocessing the sparsity pattern
    std::vector<size_t> submatrixColIndices(S.getRowIndices());
    std::vector<size_t> nnzPerCol(S.getNNZPerCol());
    std::vector<size_t> colPointers(S.getColPointers());
    std::vector<size_t> submatrixRowIndices;
    csc_matrix MM;

    size_t numCols = S.getNumCols();
    for (size_t j = 0; j < numCols; j++) {
        std::vector<size_t> submatrixRowIndicesPerCol;
    }
 
    return MM;
}