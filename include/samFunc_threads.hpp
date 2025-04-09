#pragma once
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <iomanip>
#include <ranges>

#include "householderQR.hpp"
#include "mgsQR.hpp"
#include "cscMatrix.hpp"

#include <thread>
#include <vector>
#include <algorithm>
#include <cassert>
#include <functional>

// The SAM algorithm using std::thread
void SAM_std_thread(const csc_matrix<> &source, const csc_matrix<> &target, const sparsity_pattern<> &S, csc_matrix<> &MM, const int numThreads)
{
    // Construct the map Matrix
    const size_t sparsityNumCols = S.mNumCols;
    const size_t sparsityNumRows = S.mNumRows;
    const size_t sparsityNNZ = S.mNNZ;

    // Get the sparsity pattern information
    const std::vector<ptrdiff_t> &sparsityRowIndices = S.mRowIndices;
    const std::vector<ptrdiff_t> &sparsityColPointers = S.mColPointers;

    // Get the values and row indices from the source matrix
    const std::vector<double> &sourceMatrixValues = source.mValues;
    const std::vector<ptrdiff_t> &sourceMatrixRowIndices = source.mRowIndices;
    const std::vector<ptrdiff_t> &sourceMatrixColPointers = source.mColPointers;

    // Get the values and row indices from the target matrix
    const std::vector<double> &targetMatrixValues = target.mValues;
    const std::vector<ptrdiff_t> &targetMatrixRowIndices = target.mRowIndices;
    const std::vector<ptrdiff_t> &targetMatrixColPointers = target.mColPointers;

    // Values vector for the map matrix MM
    std::vector<double> mapValues(sparsityNNZ);

    // Number of threads
    std::vector<std::thread> threads(numThreads);

    // Lambda to process a range of columns
    auto processRange = [&](ptrdiff_t colStart, ptrdiff_t colEnd)
    {
        for (ptrdiff_t j = colStart; j < colEnd; ++j)
        {
            std::vector<int> marker(sparsityNumRows, -1);
            std::vector<ptrdiff_t> J;
            std::vector<ptrdiff_t> I;
            std::vector<double> rhs;
            std::vector<std::vector<double>> submatrix;
            std::vector<double> mapColumn;

            const ptrdiff_t colBeg = sparsityColPointers[j];
            const ptrdiff_t colEnd = sparsityColPointers[j + 1];

            J.reserve(colEnd - colBeg);
            I.reserve(2 * (colEnd - colBeg));

            J.assign(sparsityRowIndices.begin() + colBeg, sparsityRowIndices.begin() + colEnd);
            assert(static_cast<ptrdiff_t>(J.size()) == (colEnd - colBeg) && "Unexpected column dimension of submatrix");

            for (ptrdiff_t i = colBeg; i < colEnd; ++i)
            {
                ptrdiff_t rowIdx = sparsityRowIndices[i];

                for (ptrdiff_t start = sparsityColPointers[rowIdx], end = sparsityColPointers[rowIdx + 1]; start < end; ++start)
                {
                    ptrdiff_t submatrixRowIdx = sparsityRowIndices[start];
                    if (marker[submatrixRowIdx] < 0)
                    {
                        marker[submatrixRowIdx] = 1;
                        I.push_back(submatrixRowIdx);
                    }
                }
            }

            std::sort(I.begin(), I.end());

            rhs.assign(I.size(), 0.0);

            ptrdiff_t targetColBeg = targetMatrixColPointers[j];
            ptrdiff_t targetColEnd = targetMatrixColPointers[j + 1];

            for (ptrdiff_t i = 0, targetRow = targetColBeg; (i < static_cast<ptrdiff_t>(I.size())) && (targetRow < targetColEnd);)
            {
                if (I[i] == targetMatrixRowIndices[targetRow])
                {
                    rhs[i] = targetMatrixValues[targetRow];
                    i++;
                    targetRow++;
                    continue;
                }
                else if (I[i] > targetMatrixRowIndices[targetRow])
                {
                    targetRow++;
                    continue;
                }
                i++;
            }

            submatrix.assign(J.size(), std::vector<double>(I.size(), 0.0));

            for (ptrdiff_t submatrixColIdx = 0, i = colBeg; i < colEnd; ++i, ++submatrixColIdx)
            {
                ptrdiff_t rowIdx = sparsityRowIndices[i];

                ptrdiff_t sourceColBeg = sourceMatrixColPointers[rowIdx];
                ptrdiff_t sourceColEnd = sourceMatrixColPointers[rowIdx + 1];

                for (ptrdiff_t k{0}, sourceRow{sourceColBeg}; (k < static_cast<ptrdiff_t>(I.size())) && (sourceRow < sourceColEnd);)
                {
                    if (I[k] == sourceMatrixRowIndices[sourceRow])
                    {
                        submatrix[submatrixColIdx][k] = sourceMatrixValues[sourceRow];
                        k++;
                        sourceRow++;
                        continue;
                    }
                    else if (I[k] > sourceMatrixRowIndices[sourceRow])
                    {
                        sourceRow++;
                        continue;
                    }
                    k++;
                }
            }

            mapColumn.assign(J.size(), 0.0);

            mgsQRSolve(submatrix, rhs, mapColumn, I.size(), J.size());

            for (ptrdiff_t i = colBeg, k = 0; i != colEnd; ++i, ++k)
            {
                mapValues[i] = mapColumn[k];
            }
        }
    };

    // Divide work among threads
    const ptrdiff_t colsPerThread = sparsityNumCols / numThreads;
    for (size_t t = 0; t < numThreads; ++t)
    {
        ptrdiff_t startCol = t * colsPerThread;
        ptrdiff_t endCol = (t == numThreads - 1) ? sparsityNumCols : (t + 1) * colsPerThread;

        threads.emplace_back(processRange, startCol, endCol);
    }

    // Join threads
    for (auto &thread : threads)
    {
        if (thread.joinable())
        {
            thread.join();
        }
    }

    MM.mNumCols = sparsityNumCols;
    MM.mNumRows = sparsityNumRows;
    MM.mNNZ = sparsityNNZ;
    MM.mRowIndices = std::ref(S.mRowIndices);
    MM.mColPointers = std::ref(S.mColPointers);
    MM.mValues = std::move(mapValues);
}