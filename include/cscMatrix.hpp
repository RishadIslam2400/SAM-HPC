#pragma once

#include <vector>
#include <iostream>
#include <numeric>

// Representation of the sparse matrices
// It uses compressed sparse column (CSC) storage format

// TODO: Refine and optimizecsc_matrix using amgcl
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

        mColPointers.reserve(ncols + 1);
        mRowIndices.reserve(mNNZ);
        mValues.reserve(mNNZ);

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

        mColPointers.reserve(ncols + 1);
        mRowIndices.reserve(mNNZ);
        mValues.reserve(mNNZ);

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

        mColPointers.reserve(ncols + 1);
        mRowIndices.reserve(mNNZ);
        mValues.reserve(mNNZ);

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

        mColPointers.reserve(ncols + 1);
        mRowIndices.reserve(mNNZ);
        mValues.reserve(mNNZ);

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

            val_type val() const {
                return *m_val;
            }

        private:
            const row_type *m_row;
            const row_type *m_end;
            const val_type *m_val;
    };

    col_iterator col_begin(size_t col) const {
        ptr_type p = mColPointers[col];
        ptr_type e = mColPointers[col + 1];
        return col_iterator(mRowIndices.data() + p, mRowIndices.data() + e, mValues.data() + p);
    }

    void printMatrix() const {
        for (size_t col = 0; col < mNumCols; ++col)
        {
            for (size_t idx = mColPointers[col]; idx < mColPointers[col + 1]; ++idx)
            {
                std::cout << "(" << mRowIndices[idx] + 1 << ", " << col + 1 << ") = " << mValues[idx] << std::endl;
            }
        }
    }
};