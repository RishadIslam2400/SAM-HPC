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

#include "SparseMatrixExceptions.hpp"

namespace SparseMatrix
{
    template <typename T>
    requires std::is_arithmetic_v<T>
    class CSRMatrix
    {
    public:
        // ============================= CONSTRUCTOR / DESTRUCTOR =============================================
        CSRMatrix() = delete;
        CSRMatrix(int n);        // square matrix nxn
        CSRMatrix(int m, int n); // general matrix mxn
        CSRMatrix(int m, int n, const std::vector<T> &vals, const std::vector<int> &row_pointers, const std::vector<int> &col_indices);

        // copy constructor
        CSRMatrix(const CSRMatrix<T> &other);
        CSRMatrix<T> &operator=(const CSRMatrix<T> &other);

        // move constructor
        CSRMatrix(CSRMatrix<T> &&other);
        CSRMatrix<T> &operator=(CSRMatrix<T> &&other);

        // destructor
        ~CSRMatrix();

        // ============================= GETTERS / SETTERS ========================================================
        int getRowCount() const;
        int getColumnCount() const;

        // Unsafe interfaces - can ruin the structure of CSR matrix
        const std::vector<T> &getValues();
        const std::vector<int> &getColIndices();
        const std::vector<int> &getRowPointers();

        // ============================= Values ==========================================================
        T get(int row, int col) const;
        CSRMatrix<T> &set(T val, int row, int col);

        // ============================= Matrix Operations ==========================================================
        std::vector<T> multiply(const std::vector<T> &x) const;
        std::vector<T> operator*(const std::vector<T> &x) const;

        CSRMatrix<T> multiply(const CSRMatrix<T> &other) const;
        CSRMatrix<T> operator*(const CSRMatrix<T> &other) const;

        CSRMatrix<T> add(const CSRMatrix<T> &other) const;
        CSRMatrix<T> operator+(const CSRMatrix<T> &other) const;

        CSRMatrix<T> subtract(const CSRMatrix<T> &other) const;
        CSRMatrix<T> operator-(const CSRMatrix<T> &other) const;

        // ============================= Friend Functions ==========================================================
        template <typename X>
        friend bool operator==(const CSRMatrix<X> &lhs, const CSRMatrix<X> &rhs);

        template <typename X>
        friend bool operator!=(const CSRMatrix<X> &lhs, const CSRMatrix<X> &rhs);

        template <typename X>
        friend std::ostream &operator<<(std::ostream &os, const CSRMatrix<X> &matrix);

    protected:
        int row_num, col_num;
        std::vector<T> *vals;
        std::vector<int> *row_pointers, *col_indices;

        // ============================= HELPERS / VALIDATORS ==============================================
        void construct(int row_num, int col_num);
        void construct(int row_num, int col_num,
                       const std::vector<T> &vals,
                       const std::vector<int> &row_pointers,
                       const std::vector<int> &col_indices);
        void destruct();
        void deepCopy(const CSRMatrix<T> &other);
        void shallowCopy(CSRMatrix<T> &&other);
        void validateCoordinates(int row, int col) const;
        void insert(int index, int row, int col, T val);
        void remove(int index, int row);
    };

    // ============================= Constructors/Assignment/Destructor ==========================================================
    template <typename T>
    CSRMatrix<T>::CSRMatrix(int n)
    {
        this->construct(n, n);
    }

    template <typename T>
    CSRMatrix<T>::CSRMatrix(int m, int n)
    {
        this->construct(m, n);
    }

    template <typename T>
    CSRMatrix<T>::CSRMatrix(int m, int n, const std::vector<T> &vals, const std::vector<int> &row_pointers, const std::vector<int> &col_indices)
    {
        this->construct(m, n, vals, row_pointers, col_indices);
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
    CSRMatrix<T>::~CSRMatrix()
    {
        this->destruct();
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

    // ============================= Setters/Getters ==========================================================
    template <typename T>
    int CSRMatrix<T>::getRowCount() const
    {
        return this->row_num;
    }

    template <typename T>
    int CSRMatrix<T>::getColumnCount() const
    {
        return this->col_num;
    }

    template <typename T>
    inline const std::vector<T> &CSRMatrix<T>::getValues()
    {
        return *vals;
    }

    template <typename T>
    inline const std::vector<int> &CSRMatrix<T>::getColIndices()
    {
        return *col_indices;
    }

    template <typename T>
    inline const std::vector<int> &CSRMatrix<T>::getRowPointers()
    {
        return *row_pointers;
    }

    template <typename T>
    T CSRMatrix<T>::get(int row, int col) const
    {
        this->validateCoordinates(row, col);

        int current_col{-1};
        for (int pos = (*(this->row_pointers))[row]; pos < (*(this->row_pointers))[row + 1]; ++pos)
        {
            current_col = (*(this->col_indices))[pos];

            if (current_col == col)
            {
                return (*(this->vals))[pos];
            }
            else if (current_col > col)
            {
                break;
            }
        }

        return T();
    }

    template <typename T>
    CSRMatrix<T> &CSRMatrix<T>::set(T val, int row, int col)
    {
        this->validateCoordinates(row, col);

        int pos = (*(this->row_pointers))[row];
        int current_col{-1};

        for (; pos < (*(this->row_pointers))[row + 1]; ++pos)
        {
            current_col = (*(this->col_indices))[pos];

            if (current_col >= col)
            {
                break;
            }
        }

        if (current_col != col)
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
    std::vector<T> CSRMatrix<T>::multiply(const std::vector<T> &x) const
    {
        // assert((this->col_num == static_cast<int>(x.size())) && "Cannot multiply: Matrix column count and vector size do not match.");

        if (this->col_num != static_cast<int>(x.size()))
        {
            throw InvalidDimensionsException("Cannot multiply: Matrix column count and vector size do not match.");
        }

        std::vector<T> result(this->row_num, T());

        if (this->vals != nullptr)
        {
            for (int i = 0; i < this->row_num; ++i)
            {
                T sum = T();
                for (int j = (*(this->row_pointers))[i]; j < (*(this->row_pointers))[i + 1]; ++j)
                {
                    sum += (*(this->vals))[j] * x[(*(this->col_indices))[j]];
                }

                result[i] = sum;
            }
        }

        return result;
    }

    template <typename T>
    std::vector<T> CSRMatrix<T>::operator*(const std::vector<T> &x) const
    {
        return this->multiply(x);
    }

    template <typename T>
    CSRMatrix<T> CSRMatrix<T>::multiply(const CSRMatrix<T> &other) const
    {
        // assert((this->col_num == other.row_num) && "Cannot multiply: Left matrix column count and right matrix row count do not match.");

        if (this->col_num != other.row_num)
        {
            throw InvalidDimensionsException("Cannot multiply: Left matrix column count and right matrix row count do not match.");
        }

        CSRMatrix<T> result(this->row_num, other.col_num);

        T a;

        for (int i = 0; i < this->row_num; ++i)
        {
            for (int j = 0; j < other.col_num; ++j)
            {
                a = T();

                for (int k = 0; k < this->col_num; ++k)
                {
                    a += this->get(i, k) * other.get(k, j);
                }

                result.set(a, i, j);
            }
        }

        return result;
    }

    template <typename T>
    CSRMatrix<T> CSRMatrix<T>::operator*(const CSRMatrix<T> &other) const
    {
        return this->multiply(other);
    }

    template <typename T>
    CSRMatrix<T> CSRMatrix<T>::add(const CSRMatrix<T> &other) const
    {
        // assert((this->row_num == other.row_num && this->col_num == other.col_num) && "Cannot add: Matrix dimensions do not match.");

        if (this->row_num != other.row_num || this->col_num != other.col_num)
        {
            throw InvalidDimensionsException("Cannot add: Matrix dimensions do not match.");
        }

        CSRMatrix<T> result(this->row_num, this->col_num);

        for (int i = 0; i < this->row_num; ++i)
        {
            for (int j = 0; j < this->col_num; ++j)
            {
                result.set(this->get(i, j) + other.get(i, j), i, j);
            }
        }

        return result;
    }

    template <typename T>
    CSRMatrix<T> CSRMatrix<T>::operator+(const CSRMatrix<T> &other) const
    {
        return this->add(other);
    }

    template <typename T>
    CSRMatrix<T> CSRMatrix<T>::subtract(const CSRMatrix<T> &other) const
    {
        // assert((this->row_num == other.row_num && this->col_num == other.col_num) && "Cannot subtract: Matrix dimensions do not match.");

        if (this->row_num != other.row_num || this->col_num != other.col_num)
        {
            throw InvalidDimensionsException("Cannot add: Matrix dimensions do not match.");
        }

        CSRMatrix<T> result(this->row_num, this->col_num);

        for (int i = 0; i < this->row_num; ++i)
        {
            for (int j = 0; j < this->col_num; ++j)
            {
                result.set(this->get(i, j) - other.get(i, j), i, j);
            }
        }

        return result;
    }

    template <typename T>
    CSRMatrix<T> CSRMatrix<T>::operator-(const CSRMatrix<T> &other) const
    {
        return this->subtract(other);
    }

    // ============================= Helpers/Validators ==========================================================
    template <typename T>
    void CSRMatrix<T>::construct(int row_num, int col_num)
    {
        // assert((row_num >= 1 && col_num >= 1) && "Matrix dimensions cannot be zero or negative.");

        if (row_num < 1 || col_num < 1)
        {
            throw InvalidDimensionsException("Matrix dimensions cannot be zero or negative.");
        }

        this->row_num = row_num;
        this->col_num = col_num;

        this->vals = nullptr;
        this->col_indices = nullptr;
        this->row_pointers = new std::vector<int>(row_num + 1, 0);
    }

    template <typename T>
    void CSRMatrix<T>::construct(int row_num, int col_num, const std::vector<T> &vals, const std::vector<int> &row_pointers, const std::vector<int> &col_indices)
    {

        // assert((row_num >= 1 && col_num >= 1) && "Matrix dimensions cannot be zero or negative.");
        // assert((row_pointers.size() == static_cast<size_t>(row_num + 1)) && "Rows pointers array does not match matrix row dimension.");
        // assert((col_indices.size() == static_cast<size_t>(row_pointers.back())) && "Column indices array does not match nonzero count.");
        // assert((vals.size() == col_indices.size()) && "Values array does not match nonzero count.");

        if (row_num < 1 || col_num < 1)
        {
            throw InvalidDimensionsException("Matrix dimensions cannot be zero or negative.");
        }

        if (row_pointers.size() != static_cast<size_t>(row_num + 1))
        {
            throw InvalidDimensionsException("Rows pointers array does not match matrix row dimension.");
        }

        if (col_indices.size() != static_cast<size_t>(row_pointers.back()))
        {
            throw InvalidDimensionsException("Column indices array does not match nonzero count.");
        }

        if (vals.size() != col_indices.size())
        {
            throw InvalidDimensionsException("Values array does not match nonzero count.");
        }

        this->row_num = row_num;
        this->col_num = col_num;

        this->vals = new std::vector<T>(vals);
        this->col_indices = new std::vector<int>(col_indices);
        this->row_pointers = new std::vector<int>(row_pointers);
    }

    template <typename T>
    void CSRMatrix<T>::deepCopy(const CSRMatrix<T> &other)
    {
        this->row_num = other.row_num;
        this->col_num = other.col_num;
        this->row_pointers = new std::vector<int>(*(other.row_pointers));

        if (other.vals != nullptr)
        {
            this->col_indices = new std::vector<int>(*(other.col_indices));
            this->vals = new std::vector<T>(*(other.vals));
        }
    }

    template <typename T>
    void CSRMatrix<T>::validateCoordinates(int row, int col) const
    {
        // assert((row >= 0 && row < this->row_num) && "Row index out of bounds.");
        // assert((col >= 0 && col < this->col_num) && "Column index out of bounds.");

        if (row < 0 || row >= this->row_num)
        {
            throw InvalidCoordinatesException("Row index out of bounds.");
        }

        if (col < 0 || col >= this->col_num)
        {
            throw InvalidCoordinatesException("Column index out of bounds.");
        }
    }

    template <typename T>
    void CSRMatrix<T>::destruct()
    {
        if (this->vals != nullptr)
        {
            delete this->vals;
            delete this->col_indices;
        }

        delete this->row_pointers;
    }

    template <typename T>
    void CSRMatrix<T>::shallowCopy(CSRMatrix<T> &&other)
    {
        this->row_num = other.row_num;
        this->col_num = other.col_num;
        this->row_pointers = other.row_pointers;
        this->vals = other.vals;
        this->col_indices = other.col_indices;

        other.row_pointers = nullptr;
        other.vals = nullptr;
        other.col_indices = nullptr;
        other.row_num = 0;
        other.col_num = 0;
    }

    template <typename T>
    void CSRMatrix<T>::insert(int index, int row, int col, T val)
    {
        if (this->vals == nullptr)
        {
            this->vals = new std::vector<T>(1, val);
            this->col_indices = new std::vector<int>(1, col);
        }
        else
        {
            this->vals->insert(this->vals->begin() + index, val);
            this->col_indices->insert(this->col_indices->begin() + index, col);
        }

        for (int i = row; i < this->row_num; ++i)
        {
            (*(this->row_pointers))[i + 1] += 1;
        }
    }

    template <typename T>
    void CSRMatrix<T>::remove(int index, int row)
    {
        this->vals->erase(this->vals->begin() + index);
        this->col_indices->erase(this->col_indices->begin() + index);

        for (int i = row; i < this->row_num; ++i)
        {
            (*(this->row_pointers))[i + 1] -= 1;
        }
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
        for (int i = 0; i < matrix.row_num; ++i)
        {
            for (int j = 0; j < matrix.col_num; ++j)
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