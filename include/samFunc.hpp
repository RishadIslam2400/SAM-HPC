#pragma once

#include <vector>
#include <mat.h>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <iomanip>
#include <ranges>

#include "householderQR.hpp"
#include "mgsQR.hpp"

// Representation of the sparse matrices
// It uses compressed sparse column (CSC) storage format

// TODO: Make templated version of the class
// TODO: Refine and optimizecsc_matrix using amgcl
class csc_matrix
{
public:
    csc_matrix() = default;
    
    // Constructor
    csc_matrix(const std::vector<double> &values, const std::vector<size_t> &rowIndices, const std::vector<size_t> &colPointers, size_t numCols, size_t numRows, size_t nnz)
        : mValues(values), mRowIndices(rowIndices), mColPointers(colPointers), mNumCols(numCols), mNumRows(numRows), mNNZ(nnz) {}

    // Copy constructor
    csc_matrix(const csc_matrix &other)
        : mValues(other.mValues), mRowIndices(other.mRowIndices), mColPointers(other.mColPointers),
          mNumCols(other.mNumCols), mNumRows(other.mNumRows), mNNZ(other.mNNZ) {}

    // Getters - the vectors are const references 
    const std::vector<double> &getValuesRef() const { return mValues; }
    const std::vector<size_t> &getRowIndicesRef() const { return mRowIndices; }
    const std::vector<size_t> &getColPointersRef() const { return mColPointers; }
    size_t getNumCols() const { return mNumCols; }
    size_t getNumRows() const { return mNumRows; }
    size_t getNNZ() const { return mNNZ; }

    // Getters - returns the copy of the vectors
    std::vector<double> getValuesCopy() const { return mValues; }
    std::vector<size_t> getRowIndicesCopy() const { return mRowIndices; }
    std::vector<size_t> getColPointersCopy() const { return mColPointers; }

    // Setters
    void setValues(const std::vector<double> &values) { mValues = values; }
    void setValues(const std::vector<double> &&values) { mValues = std::move(values); } // Move overload

    void setRowIndices(const std::vector<size_t> &rowIndices) { mRowIndices = rowIndices; }
    void setRowIndices(const std::vector<size_t> &&rowIndices) { mRowIndices = std::move(rowIndices); } // Move overload

    void setColPointers(const std::vector<size_t> &colPointers) { mColPointers = colPointers; }
    void setColPointers(const std::vector<size_t> &&colPointers) { mColPointers = std::move(colPointers); } // Move overload

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
    size_t mNumCols;                  // total number of columns
    size_t mNumRows;                  // total number of rows
    size_t mNNZ;                      // total number of non-zero elements
};

// The SAM algorithm
csc_matrix SAM(const csc_matrix& source, const csc_matrix& target, const csc_matrix& S)
{
    csc_matrix MM;

    // Construct the map Matrix
    const size_t sparsityNumCols = S.getNumCols();
    const size_t sparsityNumRows = S.getNumRows();
    const size_t sparsityNNZ = S.getNNZ();

    std::vector<double> mapValues;
    mapValues.reserve(sparsityNNZ);

    // Get the sparsity pattern information
    const std::vector<size_t> sparsityRowIndices = S.getRowIndicesRef();
    const std::vector<size_t> sparsityColPointers = S.getColPointersRef();

    // Get the values and row indices from the source matrix
    const std::vector<double> sourceMatrixValues(source.getValuesRef());
    const std::vector<size_t> sourceMatrixRowIndices(source.getRowIndicesRef());
    const std::vector<size_t> sourceMatrixColPointers(source.getColPointersRef());

    // Get the values and row indices from the target matrix
    const std::vector<double> targetMatrixValues(target.getValuesRef());
    const std::vector<size_t> targetMatrixRowIndices(target.getRowIndicesRef());
    const std::vector<size_t> targetMatrixColPointers(target.getColPointersRef());

    // Construct the submatrix information for each column
    // std::vector<size_t> marker(sparsityNumRows, -1);
    // std::vector<size_t> J;
    // std::vector<size_t> I;
    // J.reserve(sparsityNumCols);
    // I.reserve(sparsityNNZ);

    for (size_t j = 0; j < sparsityNumCols; ++j) {
        // TODO: this is for the sequential implementation
        // TODO: Use a map insted of a vector for marker
        std::vector<int> marker(sparsityNumRows, -1);
        std::vector<size_t> J;
        std::vector<size_t> I;

        size_t colBeg = sparsityColPointers[j];
        size_t colEnd = sparsityColPointers[j + 1];

        J.reserve(colEnd - colBeg);
        I.reserve(2 * (colEnd - colBeg));

        J.assign(sparsityRowIndices.begin() + colBeg, sparsityRowIndices.begin() + colEnd);
        assert(J.size() == (colEnd - colBeg) && "Unexpected column dimension of submatrix");

        // I.clear(); // TODO: not required for sequential implementation

        // Iterate over the non zeros of the current column
        for (size_t i = colBeg; i < colEnd; ++i) {
            size_t rowIdx = sparsityRowIndices[i];

            // Iterate over the non zeros of the column specified by rowIdx
            for (size_t start = sparsityColPointers[rowIdx], end = sparsityColPointers[rowIdx + 1]; start < end; ++start) {
                size_t submatrixRowIdx = sparsityRowIndices[start];
                if (marker[submatrixRowIdx] < 0) {
                    marker[submatrixRowIdx] = 1;
                    I.push_back(submatrixRowIdx);
                }
            }
        }

        // Sort the row indices for generating the submatrix
        std::sort(I.begin(), I.end());


        // Initialize the RHS vector maintaining the dimension of the submatrix
        std::vector<double> rhs(I.size(), 0.0);

        // Extract the corresponding column from the target matrix
        size_t targetColBeg = targetMatrixColPointers[j];
        size_t targetColEnd = targetMatrixColPointers[j + 1];

        // Populate the rhs vector
        for (size_t i{0}, targetRow{targetColBeg}; (i < I.size()) && (targetRow < targetColEnd);) {
            // Use the marker to map the row index of the original matrix to the submatrix
            // marker[I[i]] = static_cast<int>(i);

            size_t submatrixRowIdx = I[i];
            size_t targetRowIdx = targetMatrixRowIndices[targetRow];
            if (submatrixRowIdx == targetRowIdx) {
                rhs[i] = targetMatrixValues[targetRow];
                // TODO: update the map/marker to keep track of the index of the corresponding row
                ++i;
                ++targetRow;
            } else if(submatrixRowIdx > targetRowIdx) {
                ++targetRow;
            } else {
                ++i;
            }
        }

        // This submatrix will be in column major order
        // TODO: make the submatrix 1D vector
        std::vector<std::vector<double>> submatrix(J.size(), std::vector<double>(I.size(), 0.0));

        // Populate the submatrix
        for (size_t submatrixColIdx = 0, i = colBeg; i < colEnd; ++i, ++submatrixColIdx) {
            size_t rowIdx = sparsityRowIndices[i];

            size_t sourceColBeg = sourceMatrixColPointers[rowIdx];
            size_t sourceColEnd = sourceMatrixColPointers[rowIdx + 1];

            for (size_t k{0}, sourceRow{sourceColBeg}; (k < I.size()) && (sourceRow < sourceColBeg);) {                
                size_t submatrixRowIdx = I[k];
                size_t sourceRowIdx = sourceMatrixRowIndices[sourceRow];
                if (submatrixRowIdx == sourceRowIdx) {
                    submatrix[submatrixColIdx][k] = sourceMatrixValues[sourceRow];
                    // TODO: update the map/marker to keep track of the index of the corresponding row
                    ++i;
                    ++sourceRow;
                } else if(submatrixRowIdx > sourceRowIdx) {
                    ++sourceRow;
                } else {
                    ++i;
                }
            }
        }

        // Initialize the map vector to store the solution
        std::vector<double> mapColumn(J.size(), 0.0);

        // Solve the submatrix using QR factorization with Householder Transformations
        // Pass the arguments by value to avoid modifying the original matrices
        householderQRSolve(submatrix, rhs, mapColumn, I.size(), J.size());

        // QR solver with modified Gram Schmidt -  both are yielding same results
        // mgsQRSolve(submatrix, rhs, mapColumn, I.size(), J.size());

        // TODO: Set the values of the map matrix
        mapValues.insert(mapValues.end(), mapColumn.begin(), mapColumn.end());
    }

    MM.setNumCols(sparsityNumCols);
    MM.setNumRows(sparsityNumRows);
    MM.setNNZ(sparsityNNZ);
    MM.setRowIndices(std::ref(S.getRowIndicesRef()));
    MM.setColPointers(std::ref(S.getColPointersRef()));
    MM.setValues(std::move(mapValues));

    return MM;
}