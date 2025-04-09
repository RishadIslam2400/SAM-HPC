#pragma once

#include <vector>
// #include <mat.h>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <iomanip>
#include <ranges>

#include "householderQR.hpp"
#include "mgsQR.hpp"
#include "cscMatrix.hpp"

// The SAM algorithm
void SAM(const csc_matrix<>& source, const csc_matrix<>& target, const sparsity_pattern<>& S, csc_matrix<>& MM)
{
    // Construct the map Matrix
    const size_t sparsityNumCols = S.mNumCols;
    const size_t sparsityNumRows = S.mNumRows;
    const size_t sparsityNNZ = S.mNNZ;

    std::vector<double> mapValues(sparsityNNZ);

    // Get the sparsity pattern information
    const std::vector<ptrdiff_t>& sparsityRowIndices = S.mRowIndices;
    const std::vector<ptrdiff_t>& sparsityColPointers = S.mColPointers;

    // Get the values and row indices from the source matrix
    const std::vector<double>& sourceMatrixValues = source.mValues;
    const std::vector<ptrdiff_t>& sourceMatrixRowIndices = source.mRowIndices;
    const std::vector<ptrdiff_t>& sourceMatrixColPointers = source.mColPointers;

    // Get the values and row indices from the target matrix
    const std::vector<double>& targetMatrixValues = target.mValues;
    const std::vector<ptrdiff_t>& targetMatrixRowIndices = target.mRowIndices;
    const std::vector<ptrdiff_t>& targetMatrixColPointers = target.mColPointers;

    // Construct the submatrix information for each column
    // std::vector<size_t> marker(sparsityNumRows, -1);
    // std::vector<size_t> J;
    // std::vector<size_t> I;
    // J.reserve(sparsityNumCols);
    // I.reserve(sparsityNNZ);

    for (ptrdiff_t j = 0; j < static_cast<ptrdiff_t>(sparsityNumCols); ++j) {
        // std::cout << std::endl;
        // std::cout << "Constructing I and J for column " << j << std::endl;
        // TODO: this is for the sequential implementation
        // TODO: Use a map insted of a vector for marker
        std::vector<int> marker(sparsityNumRows, -1); // TODO: Can it be replaced by bit manipulation?
        std::vector<ptrdiff_t> J;
        std::vector<ptrdiff_t> I;

        const ptrdiff_t colBeg = sparsityColPointers[j];
        const ptrdiff_t colEnd = sparsityColPointers[j + 1];

        J.reserve(colEnd - colBeg);
        I.reserve(2 * (colEnd - colBeg));

        J.assign(sparsityRowIndices.begin() + colBeg, sparsityRowIndices.begin() + colEnd);
        assert(static_cast<ptrdiff_t>(J.size()) == (colEnd - colBeg) && "Unexpected column dimension of submatrix");

        // I.clear(); // TODO: not required for sequential implementation

        // Iterate over the non zeros of the current column
        for (ptrdiff_t i = colBeg; i < colEnd; ++i) {
            ptrdiff_t rowIdx = sparsityRowIndices[i];

            // Iterate over the non zeros of the column specified by rowIdx
            for (ptrdiff_t start = sparsityColPointers[rowIdx], end = sparsityColPointers[rowIdx + 1]; start < end; ++start) {
                ptrdiff_t submatrixRowIdx = sparsityRowIndices[start];
                if (marker[submatrixRowIdx] < 0) {
                    marker[submatrixRowIdx] = 1;
                    I.push_back(submatrixRowIdx);
                }
            }
        }

        // Sort the row indices for generating the submatrix
        std::sort(I.begin(), I.end());

        // std::cout << "Genrerated I and J for column " << j << ", size of I: " << I.size() << ", size of J: " << J.size() << std::endl;

        // Initialize the RHS vector maintaining the dimension of the submatrix
        std::vector<double> rhs(I.size(), 0.0);

        // Extract the corresponding column from the target matrix
        ptrdiff_t targetColBeg = targetMatrixColPointers[j];
        ptrdiff_t targetColEnd = targetMatrixColPointers[j + 1];

        // Populate the rhs vector
        for (ptrdiff_t i = 0, targetRow = targetColBeg; (i < static_cast<ptrdiff_t>(I.size())) && (targetRow < targetColEnd);) {
            // Use the marker to map the row index of the original matrix to the submatrix
            // marker[I[i]] = static_cast<int>(i);

            if (I[i] == targetMatrixRowIndices[targetRow]) {
                rhs[i] = targetMatrixValues[targetRow];
                // TODO: update the map/marker to keep track of the index of the corresponding row
                i++;
                targetRow++;
                continue;
            } else if(I[i] > targetMatrixRowIndices[targetRow]) {
                targetRow++;
                continue;
            }
            i++;
        }

        // std::cout << "Genrerated the rhs vector for column " << j << ", size of rhs: " << rhs.size() << std::endl;

        // This submatrix will be in column major order
        // TODO: make the submatrix 1D vector
        std::vector<std::vector<double>> submatrix(J.size(), std::vector<double>(I.size(), 0.0));

        // std::cout << "Starting to populate the submatrix for column " << j << std::endl;

        // Populate the submatrix
        for (ptrdiff_t submatrixColIdx = 0, i = colBeg; i < colEnd; ++i, ++submatrixColIdx) {
            ptrdiff_t rowIdx = sparsityRowIndices[i];

            ptrdiff_t sourceColBeg = sourceMatrixColPointers[rowIdx];
            ptrdiff_t sourceColEnd = sourceMatrixColPointers[rowIdx + 1];

            for (ptrdiff_t k{0}, sourceRow{sourceColBeg}; (k < static_cast<ptrdiff_t>(I.size())) && (sourceRow < sourceColEnd);) {
                if (I[k] == sourceMatrixRowIndices[sourceRow]) {
                    submatrix[submatrixColIdx][k] = sourceMatrixValues[sourceRow];
                    // TODO: update the map/marker to keep track of the index of the corresponding row
                    k++;
                    sourceRow++;
                    continue;
                } else if(I[k] > sourceMatrixRowIndices[sourceRow]) {
                    sourceRow++;
                    continue;
                }
                k++;
            }
        }

        // std::cout << "Populated the submatrix for column " << j << std::endl;

        // Initialize the map vector to store the solution
        std::vector<double> mapColumn(J.size(), 0.0);

        // Solve the submatrix using QR factorization with Householder Transformations
        // Pass the arguments by value to avoid modifying the original matrices
        // householderQRSolve(submatrix, rhs, mapColumn, I.size(), J.size());

        // QR solver with modified Gram Schmidt -  both are yielding same results
        mgsQRSolve(submatrix, rhs, mapColumn, I.size(), J.size());

        // mapValues.insert(mapValues.end(), mapColumn.begin(), mapColumn.end());
        for (ptrdiff_t i = colBeg, k = 0; i != colEnd; ++i, ++k) {
            mapValues[i] = mapColumn[k];
        }
    }

    MM.mNumCols = sparsityNumCols;
    MM.mNumRows = sparsityNumRows;
    MM.mNNZ = sparsityNNZ;
    // TODO: this is implcit copy. Think about using different approach
    MM.mRowIndices = std::ref(S.mRowIndices);
    MM.mColPointers = std::ref(S.mColPointers);
    MM.mValues = std::move(mapValues);
}