#include "cscMatrix.hpp"
#include <algorithm>

// This implementation uses a random guess to approximate the number of non-zeros of the resultant matrix
// so that we can finish the computation in a single pass
// TODO: find a better prediction of the number of non-zeros
void binary_spgemm_col_by_col(const sparsity_pattern<> &A, const sparsity_pattern<> &B, sparsity_pattern<> &C) {
    // Set the size of the index arrays and values in C
    C.mNumRows = A.mNumRows;
    C.mNumCols = B.mNumCols;
    C.mColPointers.reserve(C.mNumCols + 1);
    C.mRowIndices.reserve(C.mNumRows * C.mNumCols);
    C.mColPointers[0] = 0;

    // TODO: Parallelize

    // Keep track of the rows for each column of C and find out the numbner of non-zeros in each columnn
    // TODO: sequential version
    std::vector<ptrdiff_t> marker(A.mNumRows, -1);

    ptrdiff_t C_rows_nnz = 0;


    // Iterate over all the columns of B
    for (ptrdiff_t col = 0; col < B.mNumCols; ++col) {
        
        // Iterate over the non-zeros elements of the current column of B
        for (ptrdiff_t b_row = B.mColPointers[col], b_row_end = B.mColPointers[col + 1]; b_row < b_row_end; ++b_row) {
            ptrdiff_t b_row_index = B.mRowIndices[b_row];

            // Iterate over the non-zeros of the row_index column of A
            for (ptrdiff_t a_row = A.mColPointers[b_row_index], a_row_end = A.mColPointers[b_row_index + 1]; a_row < a_row_end; ++a_row) {
                ptrdiff_t a_row_index = A.mRowIndices[a_row];
                // If the row index of C for the current column is not set, set it and increase the number of non-zeros
                // if it is set, skip it
                if (marker[a_row_index] != col) {
                    marker[a_row_index] = col;
                    C.mRowIndices[C_rows_nnz] = a_row_index;
                    ++C_rows_nnz;
                }
            }
        }
        // TODO: Sequential version
        C.mColPointers[col + 1] = C_rows_nnz;

        // Sort the columns
        // TODO: implement insertion sort function
        std::sort(C.mRowIndices.begin() + C.mColPointers[col], C.mRowIndices.begin() + C.mColPointers[col + 1]);
    }
}