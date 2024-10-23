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
    // std::cout << "Simple sparsity pattern generated" << std::endl;
}

// Generates a sparsity pattern based global threshold parameter
// Sparsifies the input matrix A based on the global threshold
// Matrix entries less than the global threshold are set to 0
// For better approximation, level 2 neighbars of the sparse matrix is taken (power of 2)
// Post sparsification, for further optimization
void sparsity_pattern_global_thresh(const csc_matrix &A,const double thresh, csc_matrix &S) {
    // std::cout << "Generating sparsity pattern based on global threshold" << std::endl;
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
    std::vector<size_t> sparsityRowIndices;
    std::vector<size_t> sparsityColPointers;
    sparsityRowIndices.reserve(nnz);
    sparsityColPointers.reserve(numCols + 1);
    sparsityColPointers.push_back(0);

    // Keep track of the nnz count per column to update colPointers array
    unsigned int nnzCount{0}; 

    for (size_t j = 0; j < numCols; ++j) {
        const size_t colStart = inputColPointers[j];
        const size_t colEnd = inputColPointers[j + 1];
        for (size_t i = colStart; i < colEnd; ++i) {
            if ((std::abs(inputVals[i]) > thresh) || (inputRowIndices[i] == j)) {
                sparsityRowIndices.push_back(inputRowIndices[i]);
                ++nnzCount;
            }
        }
        sparsityColPointers.push_back(nnzCount);
    }

    const size_t sparsityNNZ = sparsityRowIndices.size();
    const size_t sparsityNumCols = numCols;
    const size_t sparsityNumRows = numRows;
    std::vector<double> sparsityValues(sparsityNNZ, 1.0);

    // Debug Print
    // std::cout << "sparsityNNZ: " << sparsityNNZ << std::endl;
    /* std::cout << "rowIndices: " << std::endl;
    for (size_t i = 0; i < sparsityNNZ; ++i) {
        std::cout << sparsityRowIndices[i] << " ";
    }
    std::cout << std::endl;
    std::cout << "colPointers: " << std::endl;
    for (size_t i = 0; i < sparsityNumCols + 1; ++i) {
        std::cout << sparsityColPointers[i] << " ";
    }
    std::cout << std::endl; */

    // TODO: Resize the vectors if too small compared to the input vector
    // TODO: Matrix matrix multiplication of S for level 2 neighbors

    // Construct the sparsity pattern matrix
    S.setNNZ(sparsityNNZ);
    S.setNumCols(sparsityNumCols);
    S.setNumRows(sparsityNumRows);
    S.setRowIndices(std::move(sparsityRowIndices));
    S.setColPointers(std::move(sparsityColPointers));
    S.setValues(std::move(sparsityValues));

    // std::cout << "Global threshold sparsity pattern generated" << std::endl;
}

// Naive implementation
// For each column, the threshold is calculated using, thresh = max(A_col) * (1-tau)
// Matrix entries above the threshold are kept in the sparsity pattern
void sparsity_pattern_col_thresh(const csc_matrix &A, const double tau, csc_matrix &S) {
    std::vector<double> inputVals{A.getValuesRef()};
    std::vector<size_t> inputRowIndices{A.getRowIndicesRef()};
    std::vector<size_t> inputColPointers{A.getColPointersRef()};
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

    // Keep track of the nnz count per column to update colPointers array
    unsigned int nnzCount{0};

    // Sparsify each column for the sparsity pattern S
    for (size_t j = 0; j < numCols; ++j) {
        const size_t colStart = inputColPointers[j];
        const size_t colEnd = inputColPointers[j + 1];

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

// Naive implementation
// For each column, lfil largest non zero entries are kept in the sparsity pattern
// TODO: Handle how to keep the diagonal in the sparsity pattern
void sparsity_pattern_lfil_thresh(const csc_matrix &A, const unsigned int lfil, csc_matrix &S) {
    std::vector<double> inputVals{A.getValuesCopy()};
    std::vector<size_t> inputRowIndices{A.getRowIndicesCopy()};
    std::vector<size_t> inputColPointers{A.getColPointersCopy()};
    const size_t numCols = A.getNumCols();
    const size_t numRows = A.getNumRows();

    // Information to construct the sparsity pattern
    std::vector<size_t> sparsityRowIndices;
    std::vector<size_t> sparsityColPointers;
    sparsityRowIndices.reserve(lfil * numCols);
    sparsityColPointers.reserve(numCols + 1);
    sparsityColPointers.push_back(0);

    // Keep track of the nnz count per column to update colPointers array
    unsigned int nnzCount{0};

    // Iterate over the columns
    for (size_t j = 0; j < numCols; ++j) {
        const size_t colStart = inputColPointers[j];
        const size_t colEnd = inputColPointers[j + 1];
        const size_t numEntriesCurrentColumn = colEnd - colStart;

        // if lfil for the current column is greater than the number of entries in the current column
        // then keep all the entries in the current column in the sparsity pattern
        if (lfil >= numEntriesCurrentColumn) {
            sparsityRowIndices.insert(sparsityRowIndices.end(), inputRowIndices.begin() + colStart, inputRowIndices.begin() + colEnd);
            nnzCount += numEntriesCurrentColumn;
            sparsityColPointers.push_back(nnzCount);
        }
        else {
            // Sort the current column in descendin order to find the lfil largest entries
            // Keep track of the row indices in the sorted column
            // TODO: May be use ranges and views to optmize
            // TODO: Separete the sorting and integrate it with the class definition
            std::vector<size_t> sortedRowIndices(inputRowIndices.begin() + colStart, inputRowIndices.begin() + colEnd);
            std::vector<double> sortedVals(inputVals.begin() + colStart, inputVals.begin() + colEnd);
            for (size_t i = 1; i < numEntriesCurrentColumn; ++i) {
                double val = sortedVals[i];
                size_t row = sortedRowIndices[i];

                // Insertion sort
                int k = i - 1;
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
            sparsityColPointers.push_back(nnzCount);
        }
    }

    // Build the rest of the sparsity pattern info
    const size_t sparsityNNZ = sparsityRowIndices.size();
    const size_t sparsityNumCols = numCols;
    const size_t sparsityNumRows = numRows;
    std::vector<double> sparsityValues(sparsityNNZ, 1.0);

    // Debug Print
    /* std::cout << "sparsityNNZ: " << sparsityNNZ << std::endl;
    std::cout << "rowIndices: " << std::endl;
    for (size_t i = 0; i < sparsityNNZ; ++i)
    {
        std::cout << sparsityRowIndices[i] << " ";
    }
    std::cout << std::endl;
    std::cout << "colPointers: " << std::endl;
    for (size_t i = 0; i < sparsityNumCols + 1; ++i)
    {
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