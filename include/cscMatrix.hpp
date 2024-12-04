#pragma once

#include <vector>
#include <iostream>
#include <numeric>
#include <memory>
#include <cassert>

//Define the macro for assert
#define ASSERTM(condition, message) \
   do { \
      assert(condition && #message); \
   } while (0)

// Representation of the sparse matrices
// It uses compressed sparse column (CSC) storage format
template <typename val_type = double, typename row_type = ptrdiff_t, typename ptr_type = row_type> 
struct csc_matrix
{
    // TODO: Use dynamic arrays instead of vectors
    std::vector<val_type> mValues;      // Non-zero elements
    std::vector<row_type> mRowIndices;  // Row indices
    std::vector<ptr_type> mColPointers; // Pointer to the start of each column in the row indices array
    size_t mNumCols;                  // total number of columns
    size_t mNumRows;                  // total number of rows
    size_t mNNZ;                      // total number of non-zero elements

    // Default constructor
    csc_matrix() = default;

    // Constructor
    // Construct the matrix from arrays, nnz not provided
    template <class PtrRange, class RowRange, class ValRange>
    csc_matrix(size_t nrows, size_t ncols, const PtrRange &colPointers, const RowRange &rowIndices, const ValRange &values) : mNumRows(nrows), mNumCols(ncols), mNNZ(0), mValues(), mRowIndices(), mColPointers() {
        static_assert(static_cast<ptrdiff_t>(nrows + 1) == std::distance(std::begin(colPointers), std::end(colPointers)), "Column Pointers has wrong size in csc constructor");

        mNNZ = colPointers[ncols];

        static_assert(static_cast<ptrdiff_t>(mNNZ) == std::distance(std::begin(rowIndices), std::end(rowIndices)), "Row Indices has wrong size in csc constructor");
        static_assert(static_cast<ptrdiff_t>(mNNZ) == std::distance(std::begin(values), std::end(values)), "Values has wrong size in csc constructor");

        mColPointers.resize(ncols + 1);
        mRowIndices.resize(mNNZ);
        mValues.resize(mNNZ);

        mColPointers[0] = colPointers[0];
        for (ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(ncols); ++i) {
            mColPointers[i + 1] = colPointers[i + 1];
            for (auto j = colPointers[i]; j < colPointers[i + 1]; ++j) {
                mRowIndices[j] = rowIndices[j];
                mValues[j] = values[j];
            }
        }
    }

    // Construct the matrix from arrays, nnz provided
    template <class PtrRange, class RowRange, class ValRange>
    csc_matrix(size_t nrows, size_t ncols, size_t nnz, const PtrRange &colPointers, const RowRange &rowIndices, const ValRange &values) : mNumRows(nrows), mNumCols(ncols), mNNZ(nnz), mValues(), mRowIndices(), mColPointers() {
        static_assert(static_cast<ptrdiff_t>(nrows + 1) == std::distance(std::begin(colPointers), std::end(colPointers)), "Column Pointers has wrong size in csc constructor");
        static_assert(static_cast<ptrdiff_t>(mNNZ) == std::distance(std::begin(rowIndices), std::end(rowIndices)), "Row Indices has wrong size in csc constructor");
        static_assert(static_cast<ptrdiff_t>(mNNZ) == std::distance(std::begin(values), std::end(values)), "Values has wrong size in csc constructor");

        mColPointers.resize(ncols + 1);
        mRowIndices.resize(mNNZ);
        mValues.resize(mNNZ);

        mColPointers[0] = colPointers[0];
        for (ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(ncols); ++i) {
            mColPointers[i + 1] = colPointers[i + 1];
            for (auto j = colPointers[i]; j < colPointers[i + 1]; ++j) {
                mRowIndices[j] = rowIndices[j];
                mValues[j] = values[j];
            }
        }
    }

    // Construct the matrix from vectors, nnz not provided
    template <class Ptr, class Row, class Val>
    csc_matrix(size_t nrows, size_t ncols, const std::vector<Ptr> &colPointers, const std::vector<Row> &rowIndices, const std::vector<Val> &values) : mNumRows(nrows), mNumCols(ncols), mNNZ(0), mValues(), mRowIndices(), mColPointers() {
        static_assert(static_cast<ptrdiff_t>(nrows + 1) == std::distance(colPointers.begin(), colPointers.end()), "Column Pointers has wrong size in csc constructor");

        mNNZ = colPointers[ncols];

        static_assert(static_cast<ptrdiff_t>(mNNZ) == std::distance(rowIndices.begin(), rowIndices.end()), "Row Indices has wrong size in csc constructor");
        static_assert(static_cast<ptrdiff_t>(mNNZ) == std::distance(values.begin(), values.end()), "Values has wrong size in csc constructor");

        mColPointers.resize(ncols + 1);
        mRowIndices.resize(mNNZ);
        mValues.resize(mNNZ);

        mColPointers = colPointers;
        mRowIndices = rowIndices;
        mValues = values;
    }

    // Construct the matrix from vectors, nnz provided
    template <class Ptr, class Row, class Val>
    csc_matrix(size_t nrows, size_t ncols, size_t nnz, const std::vector<Ptr> &colPointers, const std::vector<Row> &rowIndices, const std::vector<Val> &values) : mNumRows(nrows), mNumCols(ncols), mNNZ(nnz), mValues(), mRowIndices(), mColPointers() {
        static_assert(static_cast<ptrdiff_t>(nrows + 1) == std::distance(colPointers.begin(), colPointers.end()), "Column Pointers has wrong size in csc constructor");
        static_assert(static_cast<ptrdiff_t>(mNNZ) == std::distance(rowIndices.begin(), rowIndices.end()), "Row Indices has wrong size in csc constructor");
        static_assert(static_cast<ptrdiff_t>(mNNZ) == std::distance(values.begin(), values.end()), "Values has wrong size in csc constructor");

        mColPointers.resize(ncols + 1);
        mRowIndices.resize(mNNZ);
        mValues.resize(mNNZ);

        mColPointers = colPointers;
        mRowIndices = rowIndices;
        mValues = values;
    }

    // Copy constructor
    csc_matrix(const csc_matrix &other)
        : mNumRows(other.mNumRows), mNumCols(other.mNumCols), mNNZ(other.mNNZ), mColPointers(other.mColPointers), mRowIndices(other.mRowIndices), mValues(other.mValues) {}

    // Move constructor
    csc_matrix(csc_matrix &&other) 
        : mNumRows(other.mNumRows), mNumCols(other.mNumCols), mNNZ(other.mNNZ), mColPointers(std::move(other.mColPointers)), mRowIndices(std::move(other.mRowIndices)), mValues(std::move(other.mValues)) {
        other.mNumRows = 0;
        other.mNumCols = 0;
        other.mNNZ = 0;
    }

    const csc_matrix& operator=(const csc_matrix& other) {
        mNumRows = other.mNumRows;
        mNumCols = other.mNumCols;
        mNNZ = other.mNNZ;
        mColPointers = other.mColPointers;
        mRowIndices = other.mRowIndices;
        mValues = other.mValues;
        return *this;
    }

    const csc_matrix& operator=(const csc_matrix&& other) {
        mNumRows = other.mNumRows;
        mNumCols = other.mNumCols;
        mNNZ = other.mNNZ;
        mColPointers = std::move(other.mColPointers);
        mRowIndices = std::move(other.mRowIndices);
        mValues = std::move(other.mValues);

        other.mNumRows = 0;
        other.mNumCols = 0;
        other.mNNZ = 0;
        return *this;
    }

    class col_iterator {
        public:
            // Constructor
            col_iterator(const row_type *row, const row_type *end, const val_type *val) : m_row(row), m_end(end), m_val(val) {}

            operator bool() const {
                return m_row < m_end;
            }

            col_iterator& operator++() {
                ++m_row;
                ++m_val;
                return *this;
            }

            row_type row() const {
                return *m_row;
            }

            val_type value() const {
                return *m_val;
            }

        private:
            const row_type *m_row;
            const row_type *m_end;
            const val_type *m_val;
    };

    // Returns an iterator to the starting of a column
    col_iterator col_begin(size_t col) const {
        ptr_type p = mColPointers[col];
        ptr_type e = mColPointers[col + 1];
        return col_iterator(mRowIndices.data() + p, mRowIndices.data() + e, mValues.data() + p);
    }

    // Extract the diagonal of a csc_matrix
    std::shared_ptr<std::vector<val_type>> extractDiagonal(bool invert_sqrt = false) const {
        ASSERTM(!mColPointers.empty(), "The matrix is empty");

        // Initialize the diagonal and wrap it with a shared pointer
        auto dia = std::make_shared<std::vector<val_type>>(mNumCols);

        // retrun inverse square root of the diagonal. Replace the zero entries with 1.
        if (invert_sqrt) {
            for (ptrdiff_t j = 0; j < static_cast<ptrdiff_t>(mNumCols); ++j) {
                for (auto a = col_begin(j); a; ++a) {
                    if (a.row() == j) {
                        val_type d = std::abs(a.value());
                        d = (d == 0.0) ? 1.0 : 1.0 / std::sqrt(d); // This is an expensive operation
                        (*dia)[j] = d;
                        break;
                    }
                }
            }
        }
        // return the diagonal as is
        else {
            for (ptrdiff_t j = 0; j < static_cast<ptrdiff_t>(mNumCols); ++j) {
                for (auto a = col_begin(j); a; ++a) {
                    if (a.row() == j) {
                        val_type d = a.value();
                        (*dia)[j] = d;
                        break;
                    }
                }
            }
        }

        return dia;
    }

    // Diagonal scaling of the matrix
    std::shared_ptr<std::vector<val_type>> diagonalScale() const {
        assert(!mColPointers.empty());
        auto diag = extractDiagonal(true);
        std::vector<val_type> scaledValues = mValues;

        // TODO: Parallelize
        // Post multiplying a matrix by a diagonal matrix,
        // [a1] = d1 * [a1] (multiplying each diagonal element with the corresponding column in the matrix)
        for (ptrdiff_t j = 0; j < static_cast<ptrdiff_t>(mNumCols); ++j) {
            ptrdiff_t colStart = mColPointers[j];
            ptrdiff_t colEnd = mColPointers[j + 1];
            for (ptrdiff_t i = colStart; i < colEnd; ++i) {
                scaledValues[i] *= (*diag)[j];
            }
        }

        // Pre multiplying a matrix by a diagonal matrix,
        // [a1]^T = d1 * [a1]^T (multiplying each diagonal element with the corresponding row in the matrix)
        for (ptrdiff_t j = 0; j < static_cast<ptrdiff_t>(mNumCols); ++j) {
            ptrdiff_t colStart = mColPointers[j];
            ptrdiff_t colEnd = mColPointers[j + 1];
            for (ptrdiff_t i = colStart; i < colEnd; ++i) {
                ptrdiff_t idx = mRowIndices[i];
                scaledValues[i] *= (*diag)[idx]; // instead of iterating over the row, multiply with corresponding diagonal element
            }
        }

        auto scaledValPtr = std::make_shared<std::vector<val_type>>(scaledValues);
        return scaledValPtr;
    }

    void printMatrix() const {
        if (mColPointers.empty()) {
            std::cout << "The matrix is empty" << std::endl;
            return;
        }

        for (size_t col = 0; col < mNumCols; ++col) {
            for (size_t idx = mColPointers[col]; idx < mColPointers[col + 1]; ++idx) {
                std::cout << "(" << mRowIndices[idx] + 1 << ", " << col + 1 << ") = " << mValues[idx] << std::endl;
            }
        }
    }

    void printColumn(ptrdiff_t col_idx) const {
        if (mColPointers.empty()) {
            std::cout << "The matrix is empty" << std::endl;
            return;
        }

        ASSERTM(col_idx <= static_cast<ptrdiff_t>(mNumCols), "The column index is out of bounds");
        ptrdiff_t col_beg = mColPointers[col_idx - 1];
        ptrdiff_t col_end = mColPointers[col_idx];
        std::cout << "Column " << col_idx << ":\n";
        std::cout << "\tRow Indices: ";
        for (auto i = col_beg; i < col_end; ++i) {
            std::cout << mRowIndices[i] << " ";
        }

        std::cout << std::endl;

        std::cout << "\tValues: ";
        for (auto i = col_beg; i < col_end; ++i) {
            std::cout << mValues[i] << " ";
        }

        std::cout << std::endl << std::endl;
    }
};


// Representation of the sparsity pattern of a sparse matrix
// It uses compressed sparse column (CSC) storage format
template <typename row_type = ptrdiff_t, typename ptr_type = row_type> 
struct sparsity_pattern
{
    // TODO: Use dynamic arrays instead of vectors
    std::vector<row_type> mRowIndices;  // Row indices
    std::vector<ptr_type> mColPointers; // Pointer to the start of each column in the row indices array
    size_t mNumCols;                    // total number of columns
    size_t mNumRows;                    // total number of rows
    size_t mNNZ;                        // total number of non-zero elements

    // Default constructor
    sparsity_pattern() = default;

    // Constructor
    // Construct the matrix from arrays, nnz not provided
    template <class PtrRange, class RowRange, class ValRange>
    sparsity_pattern(size_t nrows, size_t ncols, const PtrRange &colPointers, const RowRange &rowIndices) : mNumRows(nrows), mNumCols(ncols), mNNZ(0), mRowIndices(), mColPointers() {
        static_assert(static_cast<ptrdiff_t>(nrows + 1) == std::distance(std::begin(colPointers), std::end(colPointers)), "Column Pointers has wrong size in csc constructor");

        mNNZ = colPointers[ncols];

        static_assert(static_cast<ptrdiff_t>(mNNZ) == std::distance(std::begin(rowIndices), std::end(rowIndices)), "Row Indices has wrong size in csc constructor");

        mColPointers.resize(ncols + 1);
        mRowIndices.resize(mNNZ);

        mColPointers[0] = colPointers[0];
        for (ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(ncols); ++i) {
            mColPointers[i + 1] = colPointers[i + 1];
            for (auto j = colPointers[i]; j < colPointers[i + 1]; ++j) {
                mRowIndices[j] = rowIndices[j];
            }
        }
    }

    // Construct the matrix from arrays, nnz provided
    template <class PtrRange, class RowRange, class ValRange>
    sparsity_pattern(size_t nrows, size_t ncols, size_t nnz, const PtrRange &colPointers, const RowRange &rowIndices) : mNumRows(nrows), mNumCols(ncols), mNNZ(nnz), mRowIndices(), mColPointers() {
        static_assert(static_cast<ptrdiff_t>(nrows + 1) == std::distance(std::begin(colPointers), std::end(colPointers)), "Column Pointers has wrong size in csc constructor");
        static_assert(static_cast<ptrdiff_t>(mNNZ) == std::distance(std::begin(rowIndices), std::end(rowIndices)), "Row Indices has wrong size in csc constructor");

        mColPointers.resize(ncols + 1);
        mRowIndices.resize(mNNZ);

        mColPointers[0] = colPointers[0];
        for (ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(ncols); ++i) {
            mColPointers[i + 1] = colPointers[i + 1];
            for (auto j = colPointers[i]; j < colPointers[i + 1]; ++j) {
                mRowIndices[j] = rowIndices[j];
            }
        }
    }

    // Construct the matrix from vectors, nnz not provided
    template <class Ptr, class Row, class Val>
    sparsity_pattern(size_t nrows, size_t ncols, const std::vector<Ptr> &colPointers, const std::vector<Row> &rowIndices) : mNumRows(nrows), mNumCols(ncols), mNNZ(0), mRowIndices(), mColPointers() {
        static_assert(static_cast<ptrdiff_t>(nrows + 1) == std::distance(colPointers.begin(), colPointers.end()), "Column Pointers has wrong size in csc constructor");

        mNNZ = colPointers[ncols];

        static_assert(static_cast<ptrdiff_t>(mNNZ) == std::distance(rowIndices.begin(), rowIndices.end()), "Row Indices has wrong size in csc constructor");

        mColPointers.resize(ncols + 1);
        mRowIndices.resize(mNNZ);

        mColPointers = colPointers;
        mRowIndices = rowIndices;
    }

    // Construct the matrix from vectors, nnz provided
    template <class Ptr, class Row, class Val>
    sparsity_pattern(size_t nrows, size_t ncols, size_t nnz, const std::vector<Ptr> &colPointers, const std::vector<Row> &rowIndices) : mNumRows(nrows), mNumCols(ncols), mNNZ(nnz), mRowIndices(), mColPointers() {
        static_assert(static_cast<ptrdiff_t>(nrows + 1) == std::distance(colPointers.begin(), colPointers.end()), "Column Pointers has wrong size in csc constructor");
        static_assert(static_cast<ptrdiff_t>(mNNZ) == std::distance(rowIndices.begin(), rowIndices.end()), "Row Indices has wrong size in csc constructor");

        mColPointers.resize(ncols + 1);
        mRowIndices.resize(mNNZ);

        mColPointers = colPointers;
        mRowIndices = rowIndices;
    }

    // Copy constructor
    sparsity_pattern(const sparsity_pattern &other)
        : mNumRows(other.mNumRows), mNumCols(other.mNumCols), mNNZ(other.mNNZ), mColPointers(other.mColPointers), mRowIndices(other.mRowIndices) {}

    // Move constructor
    sparsity_pattern(sparsity_pattern &&other) 
        : mNumRows(other.mNumRows), mNumCols(other.mNumCols), mNNZ(other.mNNZ), mColPointers(std::move(other.mColPointers)), mRowIndices(std::move(other.mRowIndices)) {
        other.mNumRows = 0;
        other.mNumCols = 0;
        other.mNNZ = 0;
    }

    const sparsity_pattern& operator=(const sparsity_pattern& other) {
        mNumRows = other.mNumRows;
        mNumCols = other.mNumCols;
        mNNZ = other.mNNZ;
        mColPointers = other.mColPointers;
        mRowIndices = other.mRowIndices;
        return *this;
    }

    const sparsity_pattern& operator=(const sparsity_pattern&& other) {
        mNumRows = other.mNumRows;
        mNumCols = other.mNumCols;
        mNNZ = other.mNNZ;
        mColPointers = std::move(other.mColPointers);
        mRowIndices = std::move(other.mRowIndices);
        
        other.mNumRows = 0;
        other.mNumCols = 0;
        other.mNNZ = 0;
        return *this;
    }

    class col_iterator {
        public:
            // Constructor
            col_iterator(const row_type *row, const row_type *end) : m_row(row), m_end(end) {}

            operator bool() const {
                return m_row < m_end;
            }

            col_iterator& operator++() {
                ++m_row;
                return *this;
            }

            row_type row() const {
                return *m_row;
            }

        private:
            const row_type *m_row;
            const row_type *m_end;
    };

    // Provides an iterator to the column
    col_iterator col_begin(size_t col) const {
        ptr_type p = mColPointers[col];
        ptr_type e = mColPointers[col + 1];
        return col_iterator(mRowIndices.data() + p, mRowIndices.data() + e);
    }

    // Apply binary sparse matrix matrix multiplication for the object
    // The output of the function is S = S * S.
    void get_level_2_neighbors() {
        // binary_spgemm_single_pass(this);
        binary_spgemm_double_pass(this);
    }

    void printMatrix() const {
        if (mColPointers.empty()) {
            std::cout << "The matrix is empty" << std::endl;
            return;
        }

        for (size_t col = 0; col < mNumCols; ++col) {
            for (ptrdiff_t idx = mColPointers[col]; idx < mColPointers[col + 1]; ++idx) {
                std::cout << "(" << mRowIndices[idx] + 1 << ", " << col + 1 << ")" << std::endl;
            }
        }
    }

    void printColumn(ptrdiff_t col_idx) const {
        if (mColPointers.empty()) {
            std::cout << "The matrix is empty" << std::endl;
            return;
        }

        ASSERTM(col_idx <= static_cast<ptrdiff_t>(mNumCols), "The column index is out of bounds");
        ptrdiff_t col_beg = mColPointers[col_idx - 1];
        ptrdiff_t col_end = mColPointers[col_idx];
        std::cout << "Column " << col_idx << ":\n";
        std::cout << "\tRow Indices: ";
        for (auto i = col_beg; i < col_end; ++i) {
            std::cout << mRowIndices[i] << " ";
        }

        std::cout << std::endl << std::endl;
    }
};

// This implementation uses a random guess to approximate the number of non-zeros of the resultant matrix
// so that we can finish the computation in a single pass
// TODO: find a better prediction of the number of non-zeros
void binary_spgemm_single_pass(sparsity_pattern<> *S) {
    // Set the size of the index arrays and values in C
    const size_t a_numrows = S->mNumRows;
    const size_t b_numcols = S->mNumCols;
    const size_t nnz = S->mNNZ;
    
    std::vector<ptrdiff_t> colPointers(b_numcols + 1);
    std::vector<ptrdiff_t> rowIndices;
    rowIndices.reserve(nnz * nnz); // GUESSING

    // TODO: Parallelize

    // Keep track of the rows for each column of C and find out the numbner of non-zeros in each columnn
    std::vector<ptrdiff_t> marker(a_numrows, -1);

    ptrdiff_t C_rows_nnz = 0;


    // Iterate over all the columns of B
    for (ptrdiff_t col = 0; col < static_cast<ptrdiff_t>(b_numcols); ++col) {
        
        // Iterate over the non-zeros elements of the current column of B
        for (ptrdiff_t b_row = S->mColPointers[col], b_row_end = S->mColPointers[col + 1]; b_row < b_row_end; ++b_row) {
            ptrdiff_t b_row_index = S->mRowIndices[b_row];

            // Iterate over the non-zeros of the row_index column of A
            for (ptrdiff_t a_row = S->mColPointers[b_row_index], a_row_end = S->mColPointers[b_row_index + 1]; a_row < a_row_end; ++a_row) {
                ptrdiff_t a_row_index = S->mRowIndices[a_row];
                // If the row index of C for the current column is not set, set it and increase the number of non-zeros
                // if it is set, skip it
                if (marker[a_row_index] != col) {
                    marker[a_row_index] = col;
                    rowIndices.push_back(a_row_index);
                    ++C_rows_nnz;
                }
            }
        }
        // TODO: Sequential version
        colPointers[col + 1] = C_rows_nnz;

        // Sort the columns
        // TODO: implement insertion sort function
        std::sort(rowIndices.begin() + colPointers[col], rowIndices.begin() + colPointers[col + 1]);
    }

    S->mNNZ = C_rows_nnz;
    S->mColPointers = std::move(colPointers);
    S->mRowIndices = std::move(rowIndices);
}

// This implementation computes the number of non-zeros in each column exactly and then allocates the
// the vector for storing the row indices
void binary_spgemm_double_pass(sparsity_pattern<> *S) {
    // Set the size of the index arrays and values in C
    const size_t a_numrows = S->mNumRows;
    const size_t b_numcols = S->mNumCols;
    std::vector<ptrdiff_t> colPointers(b_numcols + 1);

    // TODO: Parallelize

    // Keep track of the rows for each column of C and find out the numbner of non-zeros in each columnn
    std::vector<ptrdiff_t> marker(a_numrows, -1);

    // Iterate over all the columns of B
    for (ptrdiff_t col = 0; col < static_cast<ptrdiff_t>(b_numcols); ++col) {
        ptrdiff_t C_rows_nnz = 0;

        // Iterate over the non-zeros elements of the current column of B
        for (ptrdiff_t b_row = S->mColPointers[col], b_row_end = S->mColPointers[col + 1]; b_row < b_row_end; ++b_row) {
            ptrdiff_t b_row_index = S->mRowIndices[b_row];

            // Iterate over the non-zeros of the row_index column of A
            for (ptrdiff_t a_row = S->mColPointers[b_row_index], a_row_end = S->mColPointers[b_row_index + 1]; a_row < a_row_end; ++a_row) {
                ptrdiff_t a_row_index = S->mRowIndices[a_row];
                // If the row index of C for the current column is not set, set it and increase the number of non-zeros
                // if it is set, skip it
                if (marker[a_row_index] != col) {
                    marker[a_row_index] = col;
                    ++C_rows_nnz;
                }
            }
        }
        colPointers[col + 1] = C_rows_nnz;
    }

    // perform cumulative sum for the col pointers to find the number of non-zeros
    std::partial_sum(colPointers.cbegin(), colPointers.cend(), colPointers.begin());
    S->mNNZ = colPointers[b_numcols];
    std::vector<ptrdiff_t> rowIndices(S->mNNZ);

    marker.clear();
    marker.assign(a_numrows, -1);

    for (ptrdiff_t col = 0; col < static_cast<ptrdiff_t>(b_numcols); ++col) {
        ptrdiff_t col_beg = colPointers[col];
        ptrdiff_t col_end = col_beg;

        // Iterate over the non-zeros elements of the current column of B
        for (ptrdiff_t b_row = S->mColPointers[col], b_row_end = S->mColPointers[col + 1]; b_row < b_row_end; ++b_row) {
            ptrdiff_t b_row_index = S->mRowIndices[b_row];

            // Iterate over the non-zeros of the row_index column of A
            for (ptrdiff_t a_row = S->mColPointers[b_row_index], a_row_end = S->mColPointers[b_row_index + 1]; a_row < a_row_end; ++a_row) {
                ptrdiff_t a_row_index = S->mRowIndices[a_row];
                // If the row index of C for the current column is not set, set it and increase the number of non-zeros
                // if it is set, skip it
                if (marker[a_row_index] < col_beg) {
                    marker[a_row_index] = col_end;
                    rowIndices[col_end] = a_row_index;
                    ++col_end;
                }
            }
        }

        // Sort the columns
        // TODO: implement insertion sort function
        std::sort(rowIndices.begin() + col_beg, rowIndices.begin() + col_end);
    }

    S->mColPointers = std::move(colPointers);
    S->mRowIndices = std::move(rowIndices);
}