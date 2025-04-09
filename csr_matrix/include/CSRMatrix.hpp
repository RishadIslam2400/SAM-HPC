#pragma once
/**
 * This file is part of the SparseMatrix library
 *
 * @license  MIT
 * @author   Petr Kessler (https://kesspess.cz)
 * @link     https://github.com/uestla/Sparse-Matrix
 */

#include <vector>
#include <iostream>
#include <assert.h>
#include <concepts>
#include <type_traits>
#include <numeric>
#include <algorithm>
#include <memory>
#include <cmath>

#include "SparseMatrixExceptions.hpp"

enum class DiagonalOperation
{
    Plain,
    Abs,
    Inv,
    AbsSqrt,
    AbsInv,
    AbsInvSqrt,
};

namespace SparseMatrix
{
    template <typename T>
        requires std::is_arithmetic_v<T>
    class CSRMatrix
    {
    public:
        // ============================= Member variables ==========================================================
        size_t row_num, col_num, nnz;
        std::vector<T> *vals;
        std::vector<size_t> *row_pointers, *col_indices;

        // ============================= CONSTRUCTOR / DESTRUCTOR =============================================
        CSRMatrix() = default;
        CSRMatrix(const size_t m, const size_t n, const std::vector<T> &vals, const std::vector<size_t> &row_pointers, const std::vector<size_t> &col_indices);
        CSRMatrix(const size_t m, const size_t n, const size_t nnz, const std::vector<T> &vals, const std::vector<size_t> &row_pointers, const std::vector<size_t> &col_indices);
        CSRMatrix(const std::vector<T> &vals, const std::vector<size_t> &row_pointers, const std::vector<size_t> &col_indices);
        CSRMatrix(const size_t m, const size_t n, const T *vals, const size_t *row_pointers, const size_t *col_indices);
        CSRMatrix(const size_t m, const size_t n, const size_t nnz, const T *vals, const size_t *row_pointers, const size_t *col_indices);
        CSRMatrix(const std::vector<std::vector<T>> &matrix);

        // copy constructor
        CSRMatrix(const CSRMatrix<T> &other);
        CSRMatrix<T> &operator=(const CSRMatrix<T> &other);

        // move constructor
        CSRMatrix(CSRMatrix<T> &&other);
        CSRMatrix<T> &operator=(CSRMatrix<T> &&other);

        // destructor
        ~CSRMatrix();

        // ============================= Iterators ==========================================================
        class row_iterator
        {
        public:
            row_iterator(const size_t *col, const size_t *end, const T *val)
                : m_col(col), m_end(end), m_val(val) {}

            row_iterator &operator++()
            {
                ++m_col;
                ++m_val;
                return *this;
            }

            operator bool() const { return m_col < m_end; }
            size_t col() const { return *m_col; }
            T value() const { return *m_val; }

        private:
            const size_t *m_col;
            const size_t *m_end;
            const T *m_val;
        };
        row_iterator rowBegin(size_t row) const;

        // ============================= Values ==========================================================
        T get(size_t row, size_t col) const;
        CSRMatrix<T> &set(T val, size_t row, size_t col);

        // ============================= Matrix Operations ==========================================================
        std::shared_ptr<std::vector<T>> multiply(const std::vector<T> &x) const;
        std::shared_ptr<std::vector<T>> operator*(const std::vector<T> &x) const;

        std::shared_ptr<CSRMatrix<T>> multiply(const CSRMatrix<T> &other) const;
        std::shared_ptr<CSRMatrix<T>> operator*(const CSRMatrix<T> &other) const;

        std::shared_ptr<CSRMatrix<T>> add(const CSRMatrix<T> &other) const;
        std::shared_ptr<CSRMatrix<T>> operator+(const CSRMatrix<T> &other) const;

        std::shared_ptr<CSRMatrix<T>> subtract(const CSRMatrix<T> &other) const;
        std::shared_ptr<CSRMatrix<T>> operator-(const CSRMatrix<T> &other) const;

        // ============================= Linear Algebra Operations ===============================================
        void sortRows(); // sort rows of the matrix column-wise
        std::shared_ptr<CSRMatrix<T>> transpose() const;
        std::shared_ptr<std::vector<T>> diagonal(DiagonalOperation op = DiagonalOperation::Plain) const;

        // ============================= Friend Functions ==========================================================
        template <typename X>
        friend bool operator==(const CSRMatrix<X> &lhs, const CSRMatrix<X> &rhs);

        template <typename X>
        friend bool operator!=(const CSRMatrix<X> &lhs, const CSRMatrix<X> &rhs);

        template <typename X>
        friend std::ostream &operator<<(std::ostream &os, const CSRMatrix<X> &matrix);

        // ============================= HELPERS / VALIDATORS ==============================================
        void construct(const std::vector<T> &vals, const std::vector<size_t> &row_pointers, const std::vector<size_t> &col_indices);
        void destruct();
        void deepCopy(const CSRMatrix<T> &other);
        void shallowCopy(CSRMatrix<T> &&other);
        void validateCoordinates(size_t row, size_t col) const;
        void validateCoordinates(size_t row) const;
        void insert(size_t index, size_t row, size_t col, T val);
        void remove(size_t index, size_t row);
        size_t scanRowSize();

        // @todo: implement setSize(), setNonzeros(), scale()
    };

    // ============================= Constructors/Assignment/Destructor ==========================================================

    template <typename T>
    CSRMatrix<T>::CSRMatrix(size_t rows, size_t cols, const std::vector<T> &vals, const std::vector<size_t> &row_pointers, const std::vector<size_t> &col_indices)
    {
        // assert((rows >= 1 && cols >= 1) && "Matrix dimensions cannot be zero.");
        // assert((row_pointers.size() == rows + 1) && "Rows pointers array does not match matrix row dimension.");
        // assert((col_indices.size() == row_pointers.back()) && "Column indices array does not match nonzero count.");
        // assert((vals.size() == col_indices.size()) && "Values array does not match nonzero count.");

        if (rows < 1 || cols < 1)
            throw InvalidDimensionsException("Matrix dimensions cannot be zero.");

        if (row_pointers.size() != (rows + 1))
            throw InvalidDimensionsException("Rows pointers array does not match matrix row dimension.");

        if (col_indices.size() != row_pointers.back())
            throw InvalidDimensionsException("Column indices array does not match nonzero count.");

        if (vals.size() != col_indices.size())
            throw InvalidDimensionsException("Values array does not match nonzero count.");

        this->row_num = rows;
        this->col_num = cols;
        this->nnz = row_pointers.back();

        construct(vals, row_pointers, col_indices);
    }

    template <typename T>
    CSRMatrix<T>::CSRMatrix(const size_t rows, const size_t cols, const size_t nnz, const std::vector<T> &vals, const std::vector<size_t> &row_pointers, const std::vector<size_t> &col_indices)
    {
        // assert((rows >= 1 && cols >= 1) && "Matrix dimensions cannot be zero.");
        // assert((row_pointers.size() == rows + 1) && "Rows pointers array does not match matrix row dimension.");
        // assert((col_indices.size() == row_pointers.back()) && "Column indices array does not match nonzero count.");
        // assert((vals.size() == col_indices.size()) && "Values array does not match nonzero count.");
        // assert((nnz == row_pointers.back()) && "Nonzero count does not match rows pointers.");

        if (rows < 1 || cols < 1)
            throw InvalidDimensionsException("Matrix dimensions cannot be zero.");

        if (row_pointers.size() != (rows + 1))
            throw InvalidDimensionsException("Rows pointers array does not match matrix row dimension.");

        if (col_indices.size() != row_pointers.back())
            throw InvalidDimensionsException("Column indices array does not match nonzero count.");

        if (vals.size() != col_indices.size())
            throw InvalidDimensionsException("Values array does not match nonzero count.");

        if (nnz != row_pointers.back())
            throw InvalidDimensionsException("Nonzero count does not match rows pointers.");

        this->row_num = rows;
        this->col_num = cols;
        this->nnz = nnz;

        construct(vals, row_pointers, col_indices);
    }

    template <typename T>
    CSRMatrix<T>::CSRMatrix(const std::vector<T> &vals, const std::vector<size_t> &row_pointers, const std::vector<size_t> &col_indices)
    {
        // assert((col_indices.size() == row_pointers.back()) && "Column indices array does not match nonzero count.");
        // assert((vals.size() == col_indices.size()) && "Values array does not match nonzero count.");

        if (col_indices.size() != row_pointers.back())
            throw InvalidDimensionsException("Column indices array does not match nonzero count.");

        if (vals.size() != col_indices.size())
            throw InvalidDimensionsException("Values array does not match nonzero count.");

        this->row_num = row_pointers.size() - 1;
        this->col_num = *std::max_element(col_indices.begin(), col_indices.end()) + 1;
        this->nnz = row_pointers.back();

        construct(vals, row_pointers, col_indices);
    }

    template <typename T>
    CSRMatrix<T>::CSRMatrix(const size_t rows, const size_t cols, const T *vals, const size_t *row_pointers, const size_t *col_indices)
    {
        // assert((rows >= 1 && cols >= 1) && "Matrix dimensions cannot be zero.");
        // assert(((row_pointers != nullptr) && (col_indices != nullptr) && (vals != nullptr)) && "NULL Pointers.");

        if (rows < 1 || cols < 1)
            throw InvalidDimensionsException("Matrix dimensions cannot be zero.");

        if ((row_pointers == nullptr) || (col_indices == nullptr) || (vals == nullptr))
            throw InvalidDimensionsException("NULL Pointers.");

        this->row_num = rows;
        this->col_num = cols;
        this->nnz = row_pointers[rows];

        this->row_pointers = new std::vector<size_t>(row_pointers, row_pointers + rows + 1);
        this->col_indices = new std::vector<size_t>(col_indices, col_indices + this->nnz);
        this->vals = new std::vector<T>(vals, vals + this->nnz);
    }

    template <typename T>
    CSRMatrix<T>::CSRMatrix(const size_t rows, const size_t cols, const size_t nnz, const T *vals, const size_t *row_pointers, const size_t *col_indices)
    {
        // assert((rows >= 1 && cols >= 1) && "Matrix dimensions cannot be zero.");
        // assert(((row_pointers != nullptr) && (col_indices != nullptr) && (vals != nullptr)) && "NULL Pointers.");
        // assert((nnz == row_pointers[rows]) && "Nonzero count does not match rows pointers.");

        if (rows < 1 || cols < 1)
            throw InvalidDimensionsException("Matrix dimensions cannot be zero.");

        if ((row_pointers == nullptr) || (col_indices == nullptr) || (vals == nullptr))
            throw InvalidDimensionsException("NULL Pointers.");

        if (nnz != row_pointers[rows])
            throw InvalidDimensionsException("Nonzero count does not match rows pointers.");

        this->row_num = rows;
        this->col_num = cols;
        this->nnz = nnz;

        this->row_pointers = new std::vector<size_t>(row_pointers, row_pointers + rows + 1);
        this->col_indices = new std::vector<size_t>(col_indices, col_indices + this->nnz);
        this->vals = new std::vector<T>(vals, vals + this->nnz);
    }

    template <typename T>
    CSRMatrix<T>::CSRMatrix(const std::vector<std::vector<T>> &matrix)
    {
        this->row_num = matrix.size();
        this->col_num = matrix[0].size();
        this->row_pointers = new std::vector<size_t>(this->row_num + 1);
        (*(this->row_pointers))[0] = 0;

        // @todo: Parallelize this loop
        for (size_t i = 0; i < this->row_num; ++i)
        {
            int row_width = 0;
            for (size_t j = 0; j < this->col_num; ++j)
            {
                if (matrix[i][j] != T())
                    ++row_width;
            }
            (*(this->row_pointers))[i + 1] = row_width;
        }

        this->nnz = this->scanRowSize();
        this->vals = new std::vector<T>(this->nnz);
        this->col_indices = new std::vector<size_t>(this->nnz);

        // @todo: Parallelize this loop
        for (size_t i = 0; i < this->row_num; ++i)
        {
            size_t row_head = (*(this->row_pointers))[i];

            for (size_t j = 0; j < this->col_num; ++j)
            {
                if (matrix[i][j] != T())
                {
                    (*(this->vals))[row_head] = matrix[i][j];
                    (*(this->col_indices))[row_head] = j;
                    ++row_head;
                }
            }
        }
    }

    template <typename T>
    CSRMatrix<T>::CSRMatrix(const CSRMatrix<T> &other)
    {
        this->deepCopy(other);
    }

    template <typename T>
    CSRMatrix<T> &CSRMatrix<T>::operator=(const CSRMatrix<T> &other)
    {
        if (this != &other)
        {
            this->destruct();
            this->deepCopy(other);
        }

        return *this;
    }

    template <typename T>
    CSRMatrix<T>::CSRMatrix(CSRMatrix<T> &&other)
    {
        this->shallowCopy(std::move(other));
    }

    template <typename T>
    CSRMatrix<T> &CSRMatrix<T>::operator=(CSRMatrix<T> &&other)
    {
        if (this != &other)
        {
            this->destruct();
            this->shallowCopy(std::move(other));
        }

        return *this;
    }

    template <typename T>
    CSRMatrix<T>::~CSRMatrix()
    {
        this->destruct();
    }

    // ============================= Setters/Getters ==========================================================
    template <typename T>
    CSRMatrix<T>::row_iterator CSRMatrix<T>::rowBegin(size_t row) const
    {
        validateCoordinates(row);
        size_t start = (*(this->row_pointers))[row];
        size_t end = (*(this->row_pointers))[row + 1];
        return row_iterator(&(*(this->col_indices))[start], &(*(this->col_indices))[end], &(*(this->vals))[start]);
    }

    template <typename T>
    T CSRMatrix<T>::get(size_t row, size_t col) const
    {
        this->validateCoordinates(row, col);

        int current_col{-1};
        for (size_t pos = (*(this->row_pointers))[row]; pos < (*(this->row_pointers))[row + 1]; ++pos)
        {
            current_col = static_cast<int>((*(this->col_indices))[pos]);

            if (current_col == static_cast<int>(col))
            {
                return (*(this->vals))[pos];
            }
            else if (current_col > static_cast<int>(col))
            {
                break;
            }
        }

        return T();
    }

    template <typename T>
    CSRMatrix<T> &CSRMatrix<T>::set(T val, size_t row, size_t col)
    {
        this->validateCoordinates(row, col);

        size_t pos = (*(this->row_pointers))[row];
        int current_col{-1};

        for (; pos < (*(this->row_pointers))[row + 1]; ++pos)
        {
            current_col = static_cast<int>((*(this->col_indices))[pos]);

            if (current_col >= static_cast<int>(col))
            {
                break;
            }
        }

        if (current_col != static_cast<int>(col))
        {
            if (!(val == T()))
            {
                this->insert(pos, row, col, val);
            }
        }
        else if (val == T())
        {
            this->remove(pos, row);
        }
        else
        {
            (*(this->vals))[pos] = val;
        }

        return *this;
    }

    // ============================= Matrix Operations ==========================================================
    template <typename T>
    std::shared_ptr<std::vector<T>> CSRMatrix<T>::multiply(const std::vector<T> &x) const
    {
        // assert((this->row_pointers != nullptr && this->col_indices != nullptr && this->vals != nullptr) && "Cannot multiply: Matrix is empty.");
        // assert((this->col_num == x.size()) && "Cannot multiply: Matrix column count and vector size do not match.");

        if (this->row_pointers == nullptr || this->col_indices == nullptr || this->vals == nullptr)
            throw InvalidDimensionsException("Cannot multiply: Matrix is empty.");
        if (this->col_num != x.size())
            throw InvalidDimensionsException("Cannot multiply: Matrix column count and vector size do not match.");

        std::shared_ptr<std::vector<T>> result = std::make_shared<std::vector<T>>(this->row_num, T());

        for (size_t i = 0; i < this->row_num; ++i)
        {
            T sum = T();
            for (size_t j = (*(this->row_pointers))[i]; j < (*(this->row_pointers))[i + 1]; ++j)
            {
                sum += (*(this->vals))[j] * x[(*(this->col_indices))[j]];
            }

            (*(result))[i] = sum;
        }

        return result;
    }

    template <typename T>
    std::shared_ptr<std::vector<T>> CSRMatrix<T>::operator*(const std::vector<T> &x) const
    {
        return this->multiply(x);
    }

    template <typename T>
    std::shared_ptr<CSRMatrix<T>> CSRMatrix<T>::multiply(const CSRMatrix<T> &other) const
    {
        // assert((this->row_pointers != nullptr && this->col_indices != nullptr && this->vals != nullptr) && "Cannot multiply: Matrix is empty.");
        // assert((this->col_num == other.row_num) && "Cannot multiply: Left matrix column count and right matrix row count do not match.");

        if (this->row_pointers == nullptr || this->col_indices == nullptr || this->vals == nullptr)
            throw InvalidDimensionsException("Cannot multiply: Matrix is empty.");
        if (this->col_num != other.row_num)
            throw InvalidDimensionsException("Cannot multiply: Left matrix column count and right matrix row count do not match.");

        std::shared_ptr<CSRMatrix<T>> result = std::make_shared<CSRMatrix<T>>();

        // @todo: complete multiplication function

        return result;
    }

    template <typename T>
    std::shared_ptr<CSRMatrix<T>> CSRMatrix<T>::operator*(const CSRMatrix<T> &other) const
    {
        return this->multiply(other);
    }

    template <typename T>
    std::shared_ptr<CSRMatrix<T>> CSRMatrix<T>::add(const CSRMatrix<T> &other) const
    {
        // assert((this->row_pointers != nullptr && this->col_indices != nullptr && this->vals != nullptr) && "Cannot add: Matrix is empty.");
        // assert((this->row_num == other.row_num && this->col_num == other.col_num) && "Cannot add: Matrix dimensions do not match.");

        if (this->row_pointers == nullptr || this->col_indices == nullptr || this->vals == nullptr)
            throw InvalidDimensionsException("Cannot multiply: Matrix is empty.");
        if (this->row_num != other.row_num || this->col_num != other.col_num)
            throw InvalidDimensionsException("Cannot add: Matrix dimensions do not match.");

        std::shared_ptr<CSRMatrix<T>> result = std::make_shared<CSRMatrix<T>>();

        result->row_num = this->row_num;
        result->col_num = this->col_num;
        result->row_pointers = new std::vector<size_t>(this->row_num + 1);

        // @todo: the addition loop can be done in parallel
        // Compute non zero count per row in the result matrix
        std::vector<int> marker(result->col_num, -1);
        for (size_t i = 0; i < result->row_num; ++i)
        {
            size_t resultCols = 0;

            for (size_t j = (*(this->row_pointers))[i]; j < (*(this->row_pointers))[i + 1]; ++j)
            {
                const size_t colIdx = (*(this->col_indices))[j];
                if (marker[colIdx] != static_cast<int>(i))
                {
                    marker[colIdx] = static_cast<int>(i);
                    resultCols++;
                }
            }

            for (size_t j = (*(other.row_pointers))[i]; j < (*(other.row_pointers))[i + 1]; ++j)
            {
                const size_t colIdx = (*(other.col_indices))[j];
                if (marker[colIdx] != static_cast<int>(i))
                {
                    marker[colIdx] = static_cast<int>(i);
                    resultCols++;
                }
            }

            (*(result->row_pointers))[i + 1] = resultCols;
        }

        result->nnz =result->scanRowSize();
        result->col_indices = new std::vector<size_t>(result->nnz);
        result->vals = new std::vector<T>(result->nnz);

        // Compute the column indices and values of the result matrix
        // @todo: the addition loop can be done in parallel
        marker = std::vector<int>(result->col_num, -1);
        for (size_t i = 0; i < result->row_num; ++i)
        {
            const size_t rowBeg = (*(result->row_pointers))[i];
            size_t rowEnd = rowBeg;

            for (size_t j = (*(this->row_pointers))[i]; j < (*(this->row_pointers))[i + 1]; ++j)
            {
                size_t colIdx = (*(this->col_indices))[j];
                T value = (*(this->vals))[j];

                if (marker[colIdx] < static_cast<int>(rowBeg))
                {
                    marker[colIdx] = static_cast<int>(rowEnd);
                    (*(result->col_indices))[rowEnd] = colIdx;
                    (*(result->vals))[rowEnd] = value;
                    rowEnd++;
                }
                else
                {
                    (*(result->vals))[marker[colIdx]] += value;
                }
            }

            for (size_t j = (*(other.row_pointers))[i]; j < (*(other.row_pointers))[i + 1]; ++j)
            {
                size_t colIdx = (*(other.col_indices))[j];
                T value = (*(other.vals))[j];

                if (marker[colIdx] < static_cast<int>(rowBeg))
                {
                    marker[colIdx] = static_cast<int>(rowEnd);
                    (*(result->col_indices))[rowEnd] = colIdx;
                    (*(result->vals))[rowEnd] = value;
                    rowEnd++;
                }
                else
                {
                    (*(result->vals))[marker[colIdx]] += value;
                }
            }
        }

        result->sortRows();

        return result;
    }

    template <typename T>
    std::shared_ptr<CSRMatrix<T>> CSRMatrix<T>::operator+(const CSRMatrix<T> &other) const
    {
        return this->add(other);
    }

    template <typename T>
    std::shared_ptr<CSRMatrix<T>> CSRMatrix<T>::subtract(const CSRMatrix<T> &other) const
    {
        // assert((this->row_pointers != nullptr && this->col_indices != nullptr && this->vals != nullptr) && "Cannot subtract: Matrix is empty.");
        // assert((this->row_num == other.row_num && this->col_num == other.col_num) && "Cannot subtract: Matrix dimensions do not match.");

        if (this->row_pointers == nullptr || this->col_indices == nullptr || this->vals == nullptr)
            throw InvalidDimensionsException("Cannot multiply: Matrix is empty.");
        if (this->row_num != other.row_num || this->col_num != other.col_num)
            throw InvalidDimensionsException("Cannot add: Matrix dimensions do not match.");

        std::shared_ptr<CSRMatrix<T>> result = std::make_shared<CSRMatrix<T>>();

        result->row_num = this->row_num;
        result->col_num = this->col_num;
        result->row_pointers = new std::vector<size_t>(this->row_num + 1);

        // @todo: the addition loop can be done in parallel
        // Compute non zero count per row in the result matrix
        std::vector<int> marker(result->col_num, -1);
        for (size_t i = 0; i < result->row_num; ++i)
        {
            size_t resultCols = 0;

            for (size_t j = (*(this->row_pointers))[i]; j < (*(this->row_pointers))[i + 1]; ++j)
            {
                const size_t colIdx = (*(this->col_indices))[j];
                if (marker[colIdx] != static_cast<int>(i))
                {
                    marker[colIdx] = static_cast<int>(i);
                    resultCols++;
                }
            }

            for (size_t j = (*(other.row_pointers))[i]; j < (*(other.row_pointers))[i + 1]; ++j)
            {
                const size_t colIdx = (*(other.col_indices))[j];
                if (marker[colIdx] != static_cast<int>(i))
                {
                    marker[colIdx] = static_cast<int>(i);
                    resultCols++;
                }
            }

            (*(result->row_pointers))[i + 1] = resultCols;
        }

        result->nnz = result->scanRowSize();
        result->col_indices = new std::vector<size_t>(result->nnz);
        result->vals = new std::vector<T>(result->nnz);

        // Compute the column indices and values of the result matrix
        // @todo: the addition loop can be done in parallel
        marker = std::vector<int>(result->col_num, -1);
        for (size_t i = 0; i < result->row_num; ++i)
        {
            const size_t rowBeg = (*(result->row_pointers))[i];
            size_t rowEnd = rowBeg;

            for (size_t j = (*(this->row_pointers))[i]; j < (*(this->row_pointers))[i + 1]; ++j)
            {
                size_t colIdx = (*(this->col_indices))[j];
                T value = (*(this->vals))[j];

                if (marker[colIdx] < static_cast<int>(rowBeg))
                {
                    marker[colIdx] = static_cast<int>(rowEnd);
                    (*(result->col_indices))[rowEnd] = colIdx;
                    (*(result->vals))[rowEnd] = value;
                    rowEnd++;
                }
                else
                {
                    (*(result->vals))[marker[colIdx]] -= value;
                }
            }

            for (size_t j = (*(other.row_pointers))[i]; j < (*(other.row_pointers))[i + 1]; ++j)
            {
                size_t colIdx = (*(other.col_indices))[j];
                T value = (*(other.vals))[j];

                if (marker[colIdx] < static_cast<int>(rowBeg))
                {
                    marker[colIdx] = static_cast<int>(rowEnd);
                    (*(result->col_indices))[rowEnd] = colIdx;
                    (*(result->vals))[rowEnd] = (-value);
                    rowEnd++;
                }
                else
                {
                    (*(result->vals))[marker[colIdx]] -= value;
                }
            }
        }

        result->sortRows();

        return result;
    }

    template <typename T>
    std::shared_ptr<CSRMatrix<T>> CSRMatrix<T>::operator-(const CSRMatrix<T> &other) const
    {
        return this->subtract(other);
    }

    // ============================= Linear Algebra Operations ==========================================================
    template <typename T>
    void CSRMatrix<T>::sortRows()
    {
        assert((this->row_pointers != nullptr && this->col_indices != nullptr && this->vals != nullptr) && "Cannot sort: Matrix is empty.");
        // @todo: parallel for each row
        for (size_t i = 0; i < this->row_num; ++i)
        {
            size_t start = (*(this->row_pointers))[i];
            size_t end = (*(this->row_pointers))[i + 1];

            // insertion sort
            for (size_t j = start; j < end; ++j)
            {
                size_t col_index = (*(this->col_indices))[j];
                T val = (*(this->vals))[j];

                int i = j - 1;
                while (i >= static_cast<int>(start) && (*(this->col_indices))[i] > col_index)
                {
                    (*(this->col_indices))[i + 1] = (*(this->col_indices))[i];
                    (*(this->vals))[i + 1] = (*(this->vals))[i];
                    i--;
                }

                (*(this->col_indices))[i + 1] = col_index;
                (*(this->vals))[i + 1] = val;
            }
        }
    }

    template <typename T>
    std::shared_ptr<CSRMatrix<T>> CSRMatrix<T>::transpose() const
    {
        assert((this->row_pointers != nullptr && this->col_indices != nullptr && this->vals != nullptr) && "Cannot transpose: Matrix is empty.");
        std::shared_ptr<CSRMatrix<T>> transposedMatrix = std::make_shared<CSRMatrix<T>>();
        transposedMatrix->row_num = this->col_num;
        transposedMatrix->col_num = this->row_num;
        transposedMatrix->row_pointers = new std::vector<size_t>(this->col_num + 1);
        transposedMatrix->col_indices = new std::vector<size_t>(this->nnz);
        transposedMatrix->vals = new std::vector<T>(this->nnz);

        // Count number of non-zero elements in each column
        for (size_t n = 0; n < this->nnz; ++n)
            ++((*(transposedMatrix->row_pointers))[(*(this->col_indices))[n] + 1]);

        // Calculate prefix sum
        transposedMatrix->scanRowSize();

        for (size_t i = 0; i < this->row_num; ++i)
        {
            const size_t start = (*(this->row_pointers))[i];
            const size_t end = (*(this->row_pointers))[i + 1];
            for (size_t j = start; j < end; ++j)
            {
                size_t head = (*(transposedMatrix->row_pointers))[(*(this->col_indices))[j]]++;
                (*(transposedMatrix->col_indices))[head] = i;
                (*(transposedMatrix->vals))[head] = (*(this->vals))[j];
            }
        }

        std::rotate(transposedMatrix->row_pointers->begin(), transposedMatrix->row_pointers->begin() + this->col_num, transposedMatrix->row_pointers->begin() + this->col_num + 1);
        (*(transposedMatrix->row_pointers))[0] = 0;

        return transposedMatrix;
    }

    template <typename T>
    std::shared_ptr<std::vector<T>> CSRMatrix<T>::diagonal(DiagonalOperation op) const
    {
        assert((this->row_pointers != nullptr && this->col_indices != nullptr && this->vals != nullptr) && "Cannot extract diagonal: Matrix is empty.");
        std::shared_ptr<std::vector<T>> diag = std::make_shared<std::vector<T>>(this->row_num);

        const bool needAbsolute = (op == DiagonalOperation::Abs || op == DiagonalOperation::AbsInv || op == DiagonalOperation::AbsSqrt || op == DiagonalOperation::AbsInvSqrt);
        const bool needInvert = (op == DiagonalOperation::Inv || op == DiagonalOperation::AbsInv || op == DiagonalOperation::AbsInvSqrt);
        const bool needSqrt = (op == DiagonalOperation::AbsSqrt || op == DiagonalOperation::AbsInvSqrt);

        // @todo: this loop can be done in parallel
        for (size_t i = 0; i < this->row_num; ++i)
        {
            for (auto a = this->rowBegin(i); a; ++a)
            {
                if (a.col() == i)
                {
                    T val = a.value();

                    if (needAbsolute)
                        val = std::abs(val);

                    if (needSqrt)
                        val = std::sqrt(val);

                    if (needInvert)
                        val = (val == 0) ? 0 : 1 / val; // handling division by zero

                    (*diag)[i] = val;
                    break;
                }

                if (a.col() > i)
                    break;
            }
        }

        return diag;
    }

    // ============================= Helpers/Validators ==========================================================

    template <typename T>
    void CSRMatrix<T>::construct(const std::vector<T> &vals, const std::vector<size_t> &row_pointers, const std::vector<size_t> &col_indices)
    {
        this->row_pointers = new std::vector<size_t>(row_pointers);
        this->col_indices = new std::vector<size_t>(col_indices);
        this->vals = new std::vector<T>(vals);
    }

    template <typename T>
    void CSRMatrix<T>::deepCopy(const CSRMatrix<T> &other)
    {
        assert((other.row_pointers != nullptr && other.col_indices != nullptr && other.vals != nullptr) && "Cannot copy: Input matrix is empty.");
        this->row_num = other.row_num;
        this->col_num = other.col_num;
        this->nnz = other.nnz;
        this->row_pointers = new std::vector<size_t>(*(other.row_pointers));
        this->col_indices = new std::vector<size_t>(*(other.col_indices));
        this->vals = new std::vector<T>(*(other.vals));
    }

    template <typename T>
    void CSRMatrix<T>::validateCoordinates(size_t row, size_t col) const
    {
        // assert((row >= 0 && row < this->row_num) && "Row index out of bounds.");
        // assert((col >= 0 && col < this->col_num) && "Column index out of bounds.");

        if (row >= this->row_num)
            throw InvalidCoordinatesException("Row index out of bounds.");

        if (col >= this->col_num)
            throw InvalidCoordinatesException("Column index out of bounds.");
    }

    template <typename T>
    void CSRMatrix<T>::validateCoordinates(size_t row) const
    {
        // assert((row >= 0 && row < this->row_num) && "Row index out of bounds.");

        if (row >= this->row_num)
            throw InvalidCoordinatesException("Row index out of bounds.");
    }

    template <typename T>
    void CSRMatrix<T>::destruct()
    {
        if (this->vals != nullptr && this->row_pointers != nullptr)
        {
            delete this->vals;
            delete this->col_indices;
            delete this->row_pointers;
        }
    }

    template <typename T>
    void CSRMatrix<T>::shallowCopy(CSRMatrix<T> &&other)
    {
        assert((other->row_pointers != nullptr && other->col_indices != nullptr && other->vals != nullptr) && "Cannot move: Input matrix is empty.");
        this->row_num = other.row_num;
        this->col_num = other.col_num;
        this->nnz = other.nnz;
        this->row_pointers = other.row_pointers;
        this->vals = other.vals;
        this->col_indices = other.col_indices;

        other.row_pointers = nullptr;
        other.vals = nullptr;
        other.col_indices = nullptr;
        other.row_num = 0;
        other.col_num = 0;
        other.nnz = 0;
    }

    template <typename T>
    void CSRMatrix<T>::insert(size_t index, size_t row, size_t col, T val)
    {
        assert((this->row_pointers != nullptr && this->col_indices != nullptr && this->vals != nullptr) && "Cannot insert into an empty matrix.");
        this->vals->insert(this->vals->begin() + index, val);
        this->col_indices->insert(this->col_indices->begin() + index, col);
        this->nnz += 1;

        for (size_t i = row; i < this->row_num; ++i)
        {
            (*(this->row_pointers))[i + 1] += 1;
        }
    }

    template <typename T>
    void CSRMatrix<T>::remove(size_t index, size_t row)
    {
        assert((this->row_pointers != nullptr && this->col_indices != nullptr && this->vals != nullptr) && "Cannot remove from an empty matrix.");
        this->vals->erase(this->vals->begin() + index);
        this->col_indices->erase(this->col_indices->begin() + index);
        this->nnz -= 1;

        for (size_t i = row; i < this->row_num; ++i)
        {
            (*(this->row_pointers))[i + 1] -= 1;
        }
    }

    template <typename T>
    inline size_t CSRMatrix<T>::scanRowSize()
    {
        assert((this->row_pointers != nullptr) && "Row Pointers is null.");
        std::partial_sum(this->row_pointers->begin(), this->row_pointers->end(), this->row_pointers->begin());
        return (*(this->row_pointers))[this->row_num];
    }

    // ============================= Friend Functions ==========================================================
    template <typename T>
    bool operator==(const CSRMatrix<T> &lhs, const CSRMatrix<T> &rhs)
    {
        return (((lhs.vals == nullptr && rhs.vals == nullptr) ||
                 (lhs.vals != nullptr && rhs.vals != nullptr && *(lhs.vals) == *(rhs.vals))) &&
                ((lhs.col_indices == nullptr && rhs.col_indices == nullptr) ||
                 (lhs.col_indices != nullptr && rhs.col_indices != nullptr && *(lhs.col_indices) == *(rhs.col_indices))) &&
                (*(lhs.row_pointers) == *(rhs.row_pointers)));
    }

    template <typename T>
    bool operator!=(const CSRMatrix<T> &lhs, const CSRMatrix<T> &rhs)
    {
        return !(lhs == rhs);
    }

    template <typename T>
    std::ostream &operator<<(std::ostream &os, const CSRMatrix<T> &matrix)
    {
        assert((matrix.row_pointers != nullptr && matrix.col_indices != nullptr && matrix.vals != nullptr) && "Cannot print: Matrix is empty.");
        for (size_t i = 0; i < matrix.row_num; ++i)
        {
            for (size_t j = 0; j < matrix.col_num; ++j)
            {
                if (j != 0)
                {
                    os << " ";
                }

                os << matrix.get(i, j);
            }

            if (i < matrix.row_num - 1)
            {
                os << std::endl;
            }
        }

        return os;
    }
}