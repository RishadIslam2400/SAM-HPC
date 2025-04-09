#pragma once

#include "samFunc.hpp"

#include <memory>
#include <chrono>


// Function to generate a simple sparsity pattern
// The sparsity pattern is same as the input matrix
void simple_sparsity_pattern(const csc_matrix<> &A, sparsity_pattern<> &S) {
    // TODO: this is implcit copy. Think about using different approach
    S.mRowIndices = std::ref(A.mRowIndices);
    S.mColPointers = std::ref(A.mColPointers);
    S.mNumCols = A.mNumCols;
    S.mNumRows = A.mNumRows;
    S.mNNZ = A.mNNZ;
    // std::cout << "Number of NNZ: " << S.mNNZ << std::endl;
}

// Generates a sparsity pattern based global threshold parameter
// Sparsifies the input matrix A based on the global threshold
// Matrix entries less than the global threshold are set to 0
// For better approximation, level 2 neighbars of the sparse matrix is taken (power of 2)
// TODO: Post sparsification, for further optimization
void sparsity_pattern_global_thresh(const csc_matrix<> &A,const double thresh, sparsity_pattern<> &S) {
    // std::chrono::high_resolution_clock::time_point function_start = std::chrono::high_resolution_clock::now();
    const std::vector<ptrdiff_t>& inputRowIndices = A.mRowIndices;
    const std::vector<ptrdiff_t>& inputColPointers = A.mColPointers;
    const size_t numCols = A.mNumCols;
    const size_t numRows = A.mNumRows;
    const size_t nnz = A.mNNZ;

    // This is already scaled, values guranteed to be absolute - csc matrix interface
    // std::chrono::high_resolution_clock::time_point start_1 = std::chrono::high_resolution_clock::now();
    auto scaledValues = A.diagonalScale();
    // std::chrono::high_resolution_clock::time_point end_1 = std::chrono::high_resolution_clock::now();

    // Filter the values with the threshold with the input values array
    // Genrate the val, rowIndices and colPointers for the sparsity pattern
    // TODO: Use arrays instead of vectors for better performance
    // std::chrono::high_resolution_clock::time_point start_2 = std::chrono::high_resolution_clock::now();
    std::vector<ptrdiff_t> sparsityRowIndices; // We do not know the size beforehand as a result can not use array indexing
    std::vector<ptrdiff_t> sparsityColPointers(numCols + 1); // We know the number of elements in this vector will be same as the number of column in the original matrix
    sparsityRowIndices.reserve(nnz); // Upper bound on the size of the vector
    sparsityColPointers[0] = 0;

    // Keep track of the nnz count per column to update colPointers array
    // This is similar to a global variable. TODO: Count NNZ per column and then do reduction
    unsigned int nnzCount{0};
    for (ptrdiff_t j = 0; j < static_cast<ptrdiff_t>(numCols); ++j) {
        const ptrdiff_t colStart = inputColPointers[j];
        const ptrdiff_t colEnd = inputColPointers[j + 1];
        for (ptrdiff_t i = colStart; i < colEnd; ++i) {
            if (((*scaledValues)[i] > thresh) || (inputRowIndices[i] == j)) {
                sparsityRowIndices.push_back(inputRowIndices[i]); // TODO: Alternative ways to insert element instead of push_back
                ++nnzCount;
            }
        }
        sparsityColPointers[j + 1] = nnzCount;
    }

    const size_t sparsityNNZ = sparsityRowIndices.size();

    // Construct the sparsity pattern matrix
    S.mNNZ = sparsityNNZ;
    S.mNumCols = numCols;
    S.mNumRows = numRows;
    S.mRowIndices = std::move(sparsityRowIndices);
    S.mColPointers = std::move(sparsityColPointers);
    // std::chrono::high_resolution_clock::time_point end_2 = std::chrono::high_resolution_clock::now();

    // std::chrono::high_resolution_clock::time_point function_end = std::chrono::high_resolution_clock::now();

    // std::chrono::duration<double> time_1 = std::chrono::duration_cast<std::chrono::duration<double>>(end_1 - start_1);
    /* if (function_times != nullptr) {
        function_times[0] += time_1.count();
    }

    std::chrono::duration<double> time_2 = std::chrono::duration_cast<std::chrono::duration<double>>(end_2 - start_2);
    if (function_times != nullptr) {
        function_times[1] += time_2.count();
    }

    std::chrono::duration<double> function_time = std::chrono::duration_cast<std::chrono::duration<double>>(function_end - function_start);
    if (function_times != nullptr) {
        function_times[2] += function_time.count();
    } */

    // S.get_level_2_neighbors();

    // std::cout << "Number of NNZ: " << S.mNNZ << std::endl;
}

// Naive implementation
// For each column, the threshold is calculated using, thresh = max(A_col) * (1-tau)
// Matrix entries above the threshold are kept in the sparsity pattern
// TODO: Use a lfil parameter to fix the upper limit of non zeros in each column
void sparsity_pattern_col_thresh(const csc_matrix<> &A, const double tau, sparsity_pattern<> &S) {
    const std::vector<double>& inputVals = A.mValues;
    const std::vector<ptrdiff_t>& inputRowIndices = A.mRowIndices;
    const std::vector<ptrdiff_t>& inputColPointers = A.mColPointers;
    const size_t numCols = A.mNumCols;
    const size_t numRows = A.mNumRows;
    const size_t nnz = A.mNNZ;

    // Information to construct the sparsity pattern
    std::vector<double> sparsityValues;
    std::vector<ptrdiff_t> sparsityRowIndices;
    std::vector<ptrdiff_t> sparsityColPointers(numCols + 1);
    sparsityRowIndices.reserve(nnz); // Upper bound on the size of the vector
    sparsityColPointers[0] = 0;

    // Keep track of the nnz count per column to update colPointers array
    // TODO: Use nnz cont per column to make use of parallel calculation
    // TODO: Have scan row size implementation to create the rowPtr for the sparsity pattern
    unsigned int nnzCount{0};

    // Sparsify each column for the sparsity pattern S
    // This can be parallelized for each column
    for (ptrdiff_t j = 0; j < static_cast<ptrdiff_t>(numCols); ++j) {
        const ptrdiff_t colStart = inputColPointers[j];
        const ptrdiff_t colEnd = inputColPointers[j + 1];

        // Find the maximum absolute value in the current column
        double maxVal{0.0};
        for (ptrdiff_t i = colStart; i < colEnd; ++i) {
            if (std::abs(inputVals[i]) > maxVal) {
                maxVal = std::abs(inputVals[i]);
            }
        }

        // Find the threshold for the current column
        const double colThresh{(1.0 - tau) * maxVal};

        // Sparsify the column and construct the sparsity pattern
        for (ptrdiff_t i = colStart; i < colEnd; ++i) {
            // Keep the values greater than the threshold or the diagonal entry
            if (std::abs(inputVals[i]) > colThresh || inputRowIndices[i] == j) {
                sparsityRowIndices.push_back(inputRowIndices[i]);
                ++nnzCount;
            }
        }
        sparsityColPointers[j + 1] = nnzCount;
    }

    const size_t sparsityNNZ = sparsityRowIndices.size();
    const size_t sparsityNumCols = numCols;
    const size_t sparsityNumRows = numRows;

    // TODO: sparse matrix matrix multiplication of S for level 2 neighbors

    // Construct the sparsity pattern matrix
    S.mNNZ = sparsityNNZ;
    S.mNumCols = sparsityNumCols;
    S.mNumRows = sparsityNumRows;
    S.mRowIndices = std::move(sparsityRowIndices);
    S.mColPointers = std::move(sparsityColPointers);

    // S.get_level_2_neighbors();

    // std::cout << "Number of NNZ: " << S.mNNZ << std::endl;
}

// Naive implementation
// For each column, lfil largest non zero entries are kept in the sparsity pattern
// TODO: Handle how to keep the diagonal in the sparsity pattern
void sparsity_pattern_lfil_thresh(const csc_matrix<> &A, const unsigned int lfil, sparsity_pattern<> &S) {
    std::vector<double> inputVals = A.mValues;
    std::vector<ptrdiff_t> inputRowIndices = A.mRowIndices;
    std::vector<ptrdiff_t> inputColPointers = A.mColPointers;
    const size_t numCols = A.mNumCols;
    const size_t numRows = A.mNumRows;

    // Information to construct the sparsity pattern
    std::vector<ptrdiff_t> sparsityRowIndices;
    std::vector<ptrdiff_t> sparsityColPointers(numCols + 1);
    sparsityRowIndices.reserve(lfil * numCols); // A reasonable guess on the the size of the vector
    sparsityColPointers[0] = 0;

    // Keep track of the nnz count per column to update colPointers array
    unsigned int nnzCount{0};

    // Iterate over the columns
    for (ptrdiff_t j = 0; j < static_cast<ptrdiff_t>(numCols); ++j) {
        const ptrdiff_t colStart = inputColPointers[j];
        const ptrdiff_t colEnd = inputColPointers[j + 1];
        const ptrdiff_t numEntriesCurrentColumn = colEnd - colStart;

        // if lfil for the current column is greater than the number of entries in the current column
        // then keep all the entries in the current column in the sparsity pattern
        if (lfil >= numEntriesCurrentColumn) {
            sparsityRowIndices.insert(sparsityRowIndices.end(), inputRowIndices.begin() + colStart, inputRowIndices.begin() + colEnd);
            nnzCount += numEntriesCurrentColumn;
            sparsityColPointers[j + 1] = nnzCount;
        }
        else {
            // Sort the current column in descending order to find the lfil largest entries
            // Keep track of the row indices in the sorted column
            // TODO: May be use ranges and views to optmize
            // TODO: Separete the sorting and integrate it with the class definition
            std::vector<ptrdiff_t> sortedRowIndices(inputRowIndices.begin() + colStart, inputRowIndices.begin() + colEnd);
            std::vector<double> sortedVals(inputVals.begin() + colStart, inputVals.begin() + colEnd);

            //  Create a vector of pairs to keep track of the corresponding row indices
            std::vector<std::pair<double, ptrdiff_t>> sortedPairs;
            for (ptrdiff_t i = 1; i < numEntriesCurrentColumn; ++i) {
                double val = sortedVals[i];
                ptrdiff_t row = sortedRowIndices[i];

                // Insertion sort
                int k = static_cast<int>(i) - 1;
                while (k >= 0 && sortedVals[k] < val) {
                    sortedVals[k + 1] = sortedVals[k];
                    sortedRowIndices[k + 1] = sortedRowIndices[k];
                    --k;
                }

                sortedVals[k + 1] = val;
                sortedRowIndices[k + 1] = row;
            }

            // Construct the sparsity pattern row indices based on the sorted rows
            for (size_t i = 0; i < lfil; ++i) {
                sparsityRowIndices.push_back(sortedRowIndices[i]);
            }

            // Construct the sparsity pattern column pointers
            nnzCount += lfil;
            sparsityColPointers[j + 1] = nnzCount;
        }
    }

    // Build the rest of the sparsity pattern info
    const size_t sparsityNNZ = sparsityRowIndices.size();

    // TODO: sparse matrix matrix multiplication of S for level 2 neighbors

    // Construct the sparsity pattern matrix
    S.mNNZ = sparsityNNZ;
    S.mNumCols = numCols;
    S.mNumRows = numRows;
    S.mRowIndices = std::move(sparsityRowIndices);
    S.mColPointers = std::move(sparsityColPointers);

    // S.get_level_2_neighbors();

    // std::cout << "Number of NNZ: " << S.mNNZ << std::endl;
}