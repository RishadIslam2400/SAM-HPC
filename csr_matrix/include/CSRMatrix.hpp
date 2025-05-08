#pragma once

#include <vector>
#include <iostream>
#include <assert.h>
#include <concepts>
#include <type_traits>
#include <numeric>
#include <algorithm>
#include <memory>
#include <cmath>
#include <span>
#include <thread>
#include <functional>
#include <execution>

#include "SparseMatrixExceptions.hpp"

// GLobal thread variable for parallelization
int num_threads = 8;
constexpr bool SEQUENTIAL = false;

template <typename T>
    requires std::is_arithmetic_v<T>
class CSRMatrix {
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
    class row_iterator {
    public:
        row_iterator(const size_t *col, const size_t *end, const T *val)
            : m_col(col), m_end(end), m_val(val) {}

        row_iterator &operator++() {
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
    static void sortRow(std::span<size_t> col_indices, std::span<T> vals, int size);
    std::shared_ptr<CSRMatrix<T>> transpose() const;
    std::shared_ptr<std::vector<T>> diagonal(bool forScaling = false) const;
    void scale(T factor);

    // ============================= Friend Functions ==========================================================
    template <typename X>
    friend bool operator==(const CSRMatrix<X> &lhs, const CSRMatrix<X> &rhs);

    template <typename X>
    friend bool operator!=(const CSRMatrix<X> &lhs, const CSRMatrix<X> &rhs);

    template <typename X>
    friend std::ostream &operator<<(std::ostream &os, const CSRMatrix<X> &matrix);

    template <typename X>
    friend void swap(CSRMatrix<X> &lhs, CSRMatrix<X> &rhs);

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
    bool isEmpty() const;
};

// ============================= Constructors/Assignment/Destructor ==========================================================

template <typename T>
CSRMatrix<T>::CSRMatrix(size_t rows, size_t cols, const std::vector<T> &vals, const std::vector<size_t> &row_pointers, const std::vector<size_t> &col_indices) {
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
CSRMatrix<T>::CSRMatrix(const size_t rows, const size_t cols, const size_t nnz, const std::vector<T> &vals, const std::vector<size_t> &row_pointers, const std::vector<size_t> &col_indices) {
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
CSRMatrix<T>::CSRMatrix(const std::vector<T> &vals, const std::vector<size_t> &row_pointers, const std::vector<size_t> &col_indices) {
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
CSRMatrix<T>::CSRMatrix(const size_t rows, const size_t cols, const T *vals, const size_t *row_pointers, const size_t *col_indices) {
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
CSRMatrix<T>::CSRMatrix(const size_t rows, const size_t cols, const size_t nnz, const T *vals, const size_t *row_pointers, const size_t *col_indices) {
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
CSRMatrix<T>::CSRMatrix(const std::vector<std::vector<T>> &matrix) {
    this->row_num = matrix.size();
    this->col_num = matrix[0].size();
    this->row_pointers = new std::vector<size_t>(this->row_num + 1, 0);

    if constexpr (SEQUENTIAL) {
        for (size_t i = 0; i < this->row_num; ++i) {
            int row_width = 0;
            for (size_t j = 0; j < this->col_num; ++j) {
                if (matrix[i][j] != T())
                    ++row_width;
            }
            (*(this->row_pointers))[i + 1] = row_width;
        }
    } else {
        auto f = [&](size_t start, size_t end) {
            for (size_t i = start; i < end; ++i) {
                int row_width = 0;
                for (size_t j = 0; j < this->col_num; ++j) {
                    if (matrix[i][j] != T())
                        ++row_width;
                }
                (*(this->row_pointers))[i + 1] = row_width;
            }
        };

        launchThreads(this->row_num, f);
    }

    this->nnz = this->scanRowSize();
    this->vals = new std::vector<T>(this->nnz, T());
    this->col_indices = new std::vector<size_t>(this->nnz, 0);

    if constexpr (SEQUENTIAL) {
        for (size_t i = 0; i < this->row_num; ++i) {
            size_t row_head = (*(this->row_pointers))[i];

            for (size_t j = 0; j < this->col_num; ++j) {
                if (matrix[i][j] != T()) {
                    (*(this->vals))[row_head] = matrix[i][j];
                    (*(this->col_indices))[row_head] = j;
                    ++row_head;
                }
            }
        }
    } else {
        auto f = [&](size_t start, size_t end) {
            for (size_t i = start; i < end; ++i) {
                size_t row_head = (*(this->row_pointers))[i];

                for (size_t j = 0; j < this->col_num; ++j) {
                    if (matrix[i][j] != T()) {
                        (*(this->vals))[row_head] = matrix[i][j];
                        (*(this->col_indices))[row_head] = j;
                        ++row_head;
                    }
                }
            }
        };

        launchThreads(this->row_num, f);
    }
}

template <typename T>
CSRMatrix<T>::CSRMatrix(const CSRMatrix<T> &other) {
    this->deepCopy(other);
}

template <typename T>
CSRMatrix<T> &CSRMatrix<T>::operator=(const CSRMatrix<T> &other) {
    if (this != &other) {
        this->destruct();
        this->deepCopy(other);
    }

    return *this;
}

template <typename T>
CSRMatrix<T>::CSRMatrix(CSRMatrix<T> &&other) {
    this->shallowCopy(std::move(other));
}

template <typename T>
CSRMatrix<T> &CSRMatrix<T>::operator=(CSRMatrix<T> &&other) {
    if (this != &other) {
        this->destruct();
        this->shallowCopy(std::move(other));
    }

    return *this;
}

template <typename T>
CSRMatrix<T>::~CSRMatrix() {
    this->destruct();
}

// ============================= Setters/Getters ==========================================================
template <typename T>
CSRMatrix<T>::row_iterator CSRMatrix<T>::rowBegin(size_t row) const {
    validateCoordinates(row);
    size_t start = (*(this->row_pointers))[row];
    size_t end = (*(this->row_pointers))[row + 1];
    return row_iterator(&(*(this->col_indices))[start], &(*(this->col_indices))[end], &(*(this->vals))[start]);
}

template <typename T>
T CSRMatrix<T>::get(size_t row, size_t col) const {
    this->validateCoordinates(row, col);

    int current_col{-1};
    for (size_t pos = (*(this->row_pointers))[row]; pos < (*(this->row_pointers))[row + 1]; ++pos) {
        current_col = static_cast<int>((*(this->col_indices))[pos]);

        if (current_col == static_cast<int>(col)) {
            return (*(this->vals))[pos];
        } else if (current_col > static_cast<int>(col)) {
            break;
        }
    }

    return T();
}

template <typename T>
CSRMatrix<T> &CSRMatrix<T>::set(T val, size_t row, size_t col) {
    this->validateCoordinates(row, col);

    size_t pos = (*(this->row_pointers))[row];
    int current_col{-1};

    for (; pos < (*(this->row_pointers))[row + 1]; ++pos) {
        current_col = static_cast<int>((*(this->col_indices))[pos]);

        if (current_col >= static_cast<int>(col)) {
            break;
        }
    }

    if (current_col != static_cast<int>(col)) {
        if (!(val == T())) {
            this->insert(pos, row, col, val);
        }
    } else if (val == T()) {
        this->remove(pos, row);
    } else {
        (*(this->vals))[pos] = val;
    }

    return *this;
}

// ============================= Matrix Operations ==========================================================
template <typename T>
std::shared_ptr<std::vector<T>> CSRMatrix<T>::multiply(const std::vector<T> &x) const {
    // assert((this->row_pointers != nullptr && this->col_indices != nullptr && this->vals != nullptr) && "Cannot multiply: Matrix is empty.");
    // assert((this->col_num == x.size()) && "Cannot multiply: Matrix column count and vector size do not match.");

    if (this->row_pointers == nullptr || this->col_indices == nullptr || this->vals == nullptr)
        throw InvalidDimensionsException("Cannot multiply: Matrix is empty.");
    if (this->col_num != x.size())
        throw InvalidDimensionsException("Cannot multiply: Matrix column count and vector size do not match.");

    std::shared_ptr<std::vector<T>> result = std::make_shared<std::vector<T>>(this->row_num, T());

    if constexpr(SEQUENTIAL) {
        for (size_t i = 0; i < this->row_num; ++i) {
            T sum = T();
            for (size_t j = (*(this->row_pointers))[i]; j < (*(this->row_pointers))[i + 1]; ++j) {
                sum += (*(this->vals))[j] * x[(*(this->col_indices))[j]];
            }

            (*(result))[i] = sum;
        }
    } else {
        // @note: this is a very naive parallel implementation. look for further optimizations after more background studies
        auto f = [&](size_t start, size_t end) {
            for (size_t i = start; i < end; ++i) {
                T sum = T();
                for (size_t j = (*(this->row_pointers))[i]; j < (*(this->row_pointers))[i + 1]; ++j) {
                    sum += (*(this->vals))[j] * x[(*(this->col_indices))[j]];
                }

                (*(result))[i] = sum;
            }
        };

        launchThreads(this->row_num, f);
    }

    return result;
}

template <typename T>
std::shared_ptr<std::vector<T>> CSRMatrix<T>::operator*(const std::vector<T> &x) const
{
    return this->multiply(x);
}

template <typename T>
std::shared_ptr<CSRMatrix<T>> CSRMatrix<T>::multiply(const CSRMatrix<T> &other) const {
    // assert((this->row_pointers != nullptr && this->col_indices != nullptr && this->vals != nullptr) && "Cannot multiply: Matrix is empty.");
    // assert((this->col_num == other.row_num) && "Cannot multiply: Left matrix column count and right matrix row count do not match.");

    if (this->row_pointers == nullptr || this->col_indices == nullptr || this->vals == nullptr)
        throw InvalidDimensionsException("Cannot multiply: Matrix is empty.");
    if (this->col_num != other.row_num)
        throw InvalidDimensionsException("Cannot multiply: Left matrix column count and right matrix row count do not match.");

    std::shared_ptr<CSRMatrix<T>> result = std::make_shared<CSRMatrix<T>>();

    // todo: parallelize both functions
    if (num_threads > 16) {
        spgemm_rmerge<T>(*this, other, *result);
    }
    else {
        spgemm_saad<T>(*this, other, *result);
    }

    return result;
}

template <typename T>
std::shared_ptr<CSRMatrix<T>> CSRMatrix<T>::operator*(const CSRMatrix<T> &other) const {
    return this->multiply(other);
}

template <typename T>
std::shared_ptr<CSRMatrix<T>> CSRMatrix<T>::add(const CSRMatrix<T> &other) const {
    // assert((this->row_pointers != nullptr && this->col_indices != nullptr && this->vals != nullptr) && "Cannot add: Matrix is empty.");
    // assert((this->row_num == other.row_num && this->col_num == other.col_num) && "Cannot add: Matrix dimensions do not match.");

    if (this->row_pointers == nullptr || this->col_indices == nullptr || this->vals == nullptr)
        throw InvalidDimensionsException("Cannot multiply: Matrix is empty.");
    if (this->row_num != other.row_num || this->col_num != other.col_num)
        throw InvalidDimensionsException("Cannot add: Matrix dimensions do not match.");

    std::shared_ptr<CSRMatrix<T>> result = std::make_shared<CSRMatrix<T>>();

    result->row_num = this->row_num;
    result->col_num = this->col_num;
    result->row_pointers = new std::vector<size_t>(this->row_num + 1, 0);

    // Compute non zero count per row in the result matrix
    std::vector<int> marker(result->col_num, -1);

    if constexpr (SEQUENTIAL) {
        for (size_t i = 0; i < result->row_num; ++i) {
            size_t resultCols = 0;

            for (size_t j = (*(this->row_pointers))[i]; j < (*(this->row_pointers))[i + 1]; ++j) {
                const size_t colIdx = (*(this->col_indices))[j];
                if (marker[colIdx] != static_cast<int>(i)) {
                    marker[colIdx] = static_cast<int>(i);
                    resultCols++;
                }
            }

            for (size_t j = (*(other.row_pointers))[i]; j < (*(other.row_pointers))[i + 1]; ++j) {
                const size_t colIdx = (*(other.col_indices))[j];
                if (marker[colIdx] != static_cast<int>(i)) {
                    marker[colIdx] = static_cast<int>(i);
                    resultCols++;
                }
            }

            (*(result->row_pointers))[i + 1] = resultCols;
        }
    } else {
        auto f = [&](size_t start, size_t end) {
            for (size_t i = start; i < end; ++i) {
                size_t resultCols = 0;

                for (size_t j = (*(this->row_pointers))[i]; j < (*(this->row_pointers))[i + 1]; ++j) {
                    const size_t colIdx = (*(this->col_indices))[j];
                    if (marker[colIdx] != static_cast<int>(i)) {
                        marker[colIdx] = static_cast<int>(i);
                        resultCols++;
                    }
                }

                for (size_t j = (*(other.row_pointers))[i]; j < (*(other.row_pointers))[i + 1]; ++j) {
                    const size_t colIdx = (*(other.col_indices))[j];
                    if (marker[colIdx] != static_cast<int>(i)) {
                        marker[colIdx] = static_cast<int>(i);
                        resultCols++;
                    }
                }

                (*(result->row_pointers))[i + 1] = resultCols;
            }
        };

        launchThreads(this->row_num, f);
    }

    result->nnz = result->scanRowSize();
    result->col_indices = new std::vector<size_t>(result->nnz, 0);
    result->vals = new std::vector<T>(result->nnz, T());

    // Compute the column indices and values of the result matrix
    if constexpr (SEQUENTIAL) {
        marker = std::vector<int>(result->col_num, -1);
        
        for (size_t i = 0; i < result->row_num; ++i) {
            const size_t rowBeg = (*(result->row_pointers))[i];
            size_t rowEnd = rowBeg;

            for (size_t j = (*(this->row_pointers))[i]; j < (*(this->row_pointers))[i + 1]; ++j) {
                size_t colIdx = (*(this->col_indices))[j];
                T value = (*(this->vals))[j];

                if (marker[colIdx] < static_cast<int>(rowBeg)) {
                    marker[colIdx] = static_cast<int>(rowEnd);
                    (*(result->col_indices))[rowEnd] = colIdx;
                    (*(result->vals))[rowEnd] = value;
                    rowEnd++;
                } else {
                    (*(result->vals))[marker[colIdx]] += value;
                }
            }

            for (size_t j = (*(other.row_pointers))[i]; j < (*(other.row_pointers))[i + 1]; ++j) {
                size_t colIdx = (*(other.col_indices))[j];
                T value = (*(other.vals))[j];

                if (marker[colIdx] < static_cast<int>(rowBeg)) {
                    marker[colIdx] = static_cast<int>(rowEnd);
                    (*(result->col_indices))[rowEnd] = colIdx;
                    (*(result->vals))[rowEnd] = value;
                    rowEnd++;
                } else {
                    (*(result->vals))[marker[colIdx]] += value;
                }
            }
        }
    } else {
        auto f = [&](size_t start, size_t end) {
            marker = std::vector<int>(result->col_num, -1);
        
            for (size_t i = start; i < end; ++i) {
                const size_t rowBeg = (*(result->row_pointers))[i];
                size_t rowEnd = rowBeg;

                for (size_t j = (*(this->row_pointers))[i]; j < (*(this->row_pointers))[i + 1]; ++j) {
                    size_t colIdx = (*(this->col_indices))[j];
                    T value = (*(this->vals))[j];

                    if (marker[colIdx] < static_cast<int>(rowBeg)) {
                        marker[colIdx] = static_cast<int>(rowEnd);
                        (*(result->col_indices))[rowEnd] = colIdx;
                        (*(result->vals))[rowEnd] = value;
                        rowEnd++;
                    } else {
                        (*(result->vals))[marker[colIdx]] += value;
                    }
                }

                for (size_t j = (*(other.row_pointers))[i]; j < (*(other.row_pointers))[i + 1]; ++j) {
                    size_t colIdx = (*(other.col_indices))[j];
                    T value = (*(other.vals))[j];

                    if (marker[colIdx] < static_cast<int>(rowBeg)) {
                        marker[colIdx] = static_cast<int>(rowEnd);
                        (*(result->col_indices))[rowEnd] = colIdx;
                        (*(result->vals))[rowEnd] = value;
                        rowEnd++;
                    } else {
                        (*(result->vals))[marker[colIdx]] += value;
                    }
                }
            }
        };

        launchThreads(this->row_num, f);
    }

    result->sortRows();

    return result;
}

template <typename T>
std::shared_ptr<CSRMatrix<T>> CSRMatrix<T>::operator+(const CSRMatrix<T> &other) const {
    return this->add(other);
}

template <typename T>
std::shared_ptr<CSRMatrix<T>> CSRMatrix<T>::subtract(const CSRMatrix<T> &other) const {
    // assert((this->row_pointers != nullptr && this->col_indices != nullptr && this->vals != nullptr) && "Cannot subtract: Matrix is empty.");
    // assert((this->row_num == other.row_num && this->col_num == other.col_num) && "Cannot subtract: Matrix dimensions do not match.");

    if (this->row_pointers == nullptr || this->col_indices == nullptr || this->vals == nullptr)
        throw InvalidDimensionsException("Cannot multiply: Matrix is empty.");
    if (this->row_num != other.row_num || this->col_num != other.col_num)
        throw InvalidDimensionsException("Cannot add: Matrix dimensions do not match.");

    std::shared_ptr<CSRMatrix<T>> result = std::make_shared<CSRMatrix<T>>();

    result->row_num = this->row_num;
    result->col_num = this->col_num;
    result->row_pointers = new std::vector<size_t>(this->row_num + 1, 0);

    // Compute non zero count per row in the result matrix
    if constexpr (SEQUENTIAL) {
        std::vector<int> marker(result->col_num, -1);
        for (size_t i = 0; i < result->row_num; ++i) {
            size_t resultCols = 0;

            for (size_t j = (*(this->row_pointers))[i]; j < (*(this->row_pointers))[i + 1]; ++j) {
                const size_t colIdx = (*(this->col_indices))[j];
                if (marker[colIdx] != static_cast<int>(i)) {
                    marker[colIdx] = static_cast<int>(i);
                    resultCols++;
                }
            }

            for (size_t j = (*(other.row_pointers))[i]; j < (*(other.row_pointers))[i + 1]; ++j) {
                const size_t colIdx = (*(other.col_indices))[j];
                if (marker[colIdx] != static_cast<int>(i)) {
                    marker[colIdx] = static_cast<int>(i);
                    resultCols++;
                }
            }

            (*(result->row_pointers))[i + 1] = resultCols;
        }
    } else {
        auto f = [&](size_t start, size_t end) {
            std::vector<int> marker(result->col_num, -1);
            for (size_t i = start; i < end; ++i) {
                size_t resultCols = 0;

                for (size_t j = (*(this->row_pointers))[i]; j < (*(this->row_pointers))[i + 1]; ++j) {
                    const size_t colIdx = (*(this->col_indices))[j];
                    if (marker[colIdx] != static_cast<int>(i)) {
                        marker[colIdx] = static_cast<int>(i);
                        resultCols++;
                    }
                }

                for (size_t j = (*(other.row_pointers))[i]; j < (*(other.row_pointers))[i + 1]; ++j) {
                    const size_t colIdx = (*(other.col_indices))[j];
                    if (marker[colIdx] != static_cast<int>(i)) {
                        marker[colIdx] = static_cast<int>(i);
                        resultCols++;
                    }
                }

                (*(result->row_pointers))[i + 1] = resultCols;
            }
        };

        launchThreads(this->row_num, f);
    }

    result->nnz = result->scanRowSize();
    result->col_indices = new std::vector<size_t>(result->nnz, 0);
    result->vals = new std::vector<T>(result->nnz, T());

    // Compute the column indices and values of the result matrix
    if constexpr (SEQUENTIAL) {
        std::vector<int> marker(result->col_num, -1);
        for (size_t i = 0; i < result->row_num; ++i) {
            const size_t rowBeg = (*(result->row_pointers))[i];
            size_t rowEnd = rowBeg;

            for (size_t j = (*(this->row_pointers))[i]; j < (*(this->row_pointers))[i + 1]; ++j) {
                size_t colIdx = (*(this->col_indices))[j];
                T value = (*(this->vals))[j];

                if (marker[colIdx] < static_cast<int>(rowBeg)) {
                    marker[colIdx] = static_cast<int>(rowEnd);
                    (*(result->col_indices))[rowEnd] = colIdx;
                    (*(result->vals))[rowEnd] = value;
                    rowEnd++;
                }
                else {
                    (*(result->vals))[marker[colIdx]] -= value;
                }
            }

            for (size_t j = (*(other.row_pointers))[i]; j < (*(other.row_pointers))[i + 1]; ++j) {
                size_t colIdx = (*(other.col_indices))[j];
                T value = (*(other.vals))[j];

                if (marker[colIdx] < static_cast<int>(rowBeg)) {
                    marker[colIdx] = static_cast<int>(rowEnd);
                    (*(result->col_indices))[rowEnd] = colIdx;
                    (*(result->vals))[rowEnd] = (-value);
                    rowEnd++;
                } else {
                    (*(result->vals))[marker[colIdx]] -= value;
                }
            }
        }
    } else {
        auto f = [&](size_t start, size_t end) {
            std::vector<int> marker(result->col_num, -1);
            for (size_t i = start; i < end; ++i) {
                const size_t rowBeg = (*(result->row_pointers))[i];
                size_t rowEnd = rowBeg;

                for (size_t j = (*(this->row_pointers))[i]; j < (*(this->row_pointers))[i + 1]; ++j) {
                    size_t colIdx = (*(this->col_indices))[j];
                    T value = (*(this->vals))[j];

                    if (marker[colIdx] < static_cast<int>(rowBeg)) {
                        marker[colIdx] = static_cast<int>(rowEnd);
                        (*(result->col_indices))[rowEnd] = colIdx;
                        (*(result->vals))[rowEnd] = value;
                        rowEnd++;
                    }
                    else {
                        (*(result->vals))[marker[colIdx]] -= value;
                    }
                }

                for (size_t j = (*(other.row_pointers))[i]; j < (*(other.row_pointers))[i + 1]; ++j) {
                    size_t colIdx = (*(other.col_indices))[j];
                    T value = (*(other.vals))[j];

                    if (marker[colIdx] < static_cast<int>(rowBeg)) {
                        marker[colIdx] = static_cast<int>(rowEnd);
                        (*(result->col_indices))[rowEnd] = colIdx;
                        (*(result->vals))[rowEnd] = (-value);
                        rowEnd++;
                    } else {
                        (*(result->vals))[marker[colIdx]] -= value;
                    }
                }
            }
        };

        launchThreads(this->row_num, f);
    }

    result->sortRows();

    return result;
}

template <typename T>
std::shared_ptr<CSRMatrix<T>> CSRMatrix<T>::operator-(const CSRMatrix<T> &other) const {
    return this->subtract(other);
}

// ============================= Linear Algebra Operations ==========================================================
template <typename T>
void CSRMatrix<T>::sortRows() {
    assert((this->row_pointers != nullptr && this->col_indices != nullptr && this->vals != nullptr) && "Cannot sort: Matrix is empty.");
    
    // todo: remove span, use pointers
    if constexpr (SEQUENTIAL) {
        for (size_t i = 0; i < this->row_num; ++i) {
            size_t start = (*(this->row_pointers))[i];
            size_t end = (*(this->row_pointers))[i + 1];

            std::span<size_t> col_indices_span(this->col_indices->begin() + start, end - start);
            std::span<T> vals_span(this->vals->begin() + start, end - start);
            sortRow(col_indices_span, vals_span, end - start);
        }
    } else {
        auto f = [&](size_t start, size_t end) {
            for (size_t i = start; i < end; ++i) {
                size_t row_start = (*(this->row_pointers))[i];
                size_t row_end = (*(this->row_pointers))[i + 1];

                std::span<size_t> col_indices_span(this->col_indices->begin() + row_start, row_end - row_start);
                std::span<T> vals_span(this->vals->begin() + row_start, row_end - row_start);
                sortRow(col_indices_span, vals_span, row_end - row_start);
            }
        };

        launchThreads(this->row_num, f);
    }
}

template <typename T>
void CSRMatrix<T>::sortRow(std::span<size_t> col_indices, std::span<T> vals, int size) {
    // insertion sort
    for (int j = 1; j < size; ++j) {
        size_t colIdx = col_indices[j];
        T val = vals[j];

        int i = j - 1;

        while (i >= 0 && col_indices[i] > colIdx) {
            col_indices[i + 1] = col_indices[i];
            vals[i + 1] = vals[i];
            i--;
        }

        col_indices[i + 1] = colIdx;
        vals[i + 1] = val;
    }
}

template <typename T>
std::shared_ptr<CSRMatrix<T>> CSRMatrix<T>::transpose() const {
    assert((this->row_pointers != nullptr && this->col_indices != nullptr && this->vals != nullptr) && "Cannot transpose: Matrix is empty.");
    std::shared_ptr<CSRMatrix<T>> transposedMatrix = std::make_shared<CSRMatrix<T>>();
    transposedMatrix->row_num = this->col_num;
    transposedMatrix->col_num = this->row_num;
    transposedMatrix->row_pointers = new std::vector<size_t>(this->col_num + 1, 0);
    transposedMatrix->col_indices = new std::vector<size_t>(this->nnz, 0);
    transposedMatrix->vals = new std::vector<T>(this->nnz, T());

    // Count number of non-zero elements in each column
    for (size_t n = 0; n < this->nnz; ++n) {
        ++((*(transposedMatrix->row_pointers))[(*(this->col_indices))[n] + 1]);
    }
    transposedMatrix->scanRowSize();

    if constexpr (SEQUENTIAL) {
        for (size_t i = 0; i < this->row_num; ++i) {
            const size_t start = (*(this->row_pointers))[i];
            const size_t end = (*(this->row_pointers))[i + 1];
            for (size_t j = start; j < end; ++j) {
                size_t head = (*(transposedMatrix->row_pointers))[(*(this->col_indices))[j]]++;
                (*(transposedMatrix->col_indices))[head] = i;
                (*(transposedMatrix->vals))[head] = (*(this->vals))[j];
            }
        }
    } else {
        auto f = [&](size_t start, size_t end) {
            for (size_t i = start; i < end; ++i) {
                const size_t row_start = (*(this->row_pointers))[i];
                const size_t row_end = (*(this->row_pointers))[i + 1];
                for (size_t j = row_start; j < row_end; ++j) {
                    size_t head = (*(transposedMatrix->row_pointers))[(*(this->col_indices))[j]]++;
                    (*(transposedMatrix->col_indices))[head] = i;
                    (*(transposedMatrix->vals))[head] = (*(this->vals))[j];
                }
            }
        };

        launchThreads(this->row_num, f);
    }

    std::rotate(transposedMatrix->row_pointers->begin(), transposedMatrix->row_pointers->begin() + this->col_num, transposedMatrix->row_pointers->begin() + this->col_num + 1);
    (*(transposedMatrix->row_pointers))[0] = 0;

    return transposedMatrix;
}

template <typename T>
std::shared_ptr<std::vector<T>> CSRMatrix<T>::diagonal(bool forScaling) const {
    assert((this->row_pointers != nullptr && this->col_indices != nullptr && this->vals != nullptr) && "Cannot extract diagonal: Matrix is empty.");

    if constexpr (SEQUENTIAL) {
        if (!forScaling) {
            std::shared_ptr<std::vector<T>> diag = std::make_shared<std::vector<T>>(this->row_num);
            for (size_t i = 0; i < this->row_num; ++i) {
                for (auto a = this->rowBegin(i); a; ++a) {
                    if (a.col() == i) {
                        (*diag)[i] = a.value();
                        break;
                    }
                }
            }
            return diag;
        } else {
            std::shared_ptr<std::vector<T>> diag = std::make_shared<std::vector<T>>(this->row_num, 1);
            for (size_t i = 0; i < this->row_num; ++i) {
                for (auto a = this->rowBegin(i); a; ++a) {
                    if (a.col() == i) {
                        if (a.value() != T()) {
                            (*diag)[i] = 1 / std::sqrt(std::abs(a.value()));
                            break;
                        }
                    }
                }
            }
            return diag;
        }
    } else {
        if (!forScaling) {
            std::shared_ptr<std::vector<T>> diag = std::make_shared<std::vector<T>>(this->row_num);
            auto f = [&](size_t start, size_t end) {
                for (size_t i = start; i < end; ++i) {
                    for (auto a = this->rowBegin(i); a; ++a) {
                        if (a.col() == i) {
                            (*diag)[i] = a.value();
                            break;
                        }
                    }
                }
            };

            launchThreads(this->row_num, f);
            return diag;
        } else {
            std::shared_ptr<std::vector<T>> diag = std::make_shared<std::vector<T>>(this->row_num, 1);
            auto f = [&](size_t start, size_t end) {
                for (size_t i = start; i < end; ++i) {
                    for (auto a = this->rowBegin(i); a; ++a) {
                        if (a.col() == i) {
                            if (a.value() != T()) {
                                (*diag)[i] = 1 / std::sqrt(std::abs(a.value()));
                                break;
                            }
                        }
                    }
                }
            };

            launchThreads(this->row_num, f);
            return diag;
        }
    }

    return std::make_shared<std::vector<T>>(this->row_num, T()); // this should never be executed
}

template <typename T>
void CSRMatrix<T>::scale(T factor) {
    assert((this->row_pointers != nullptr && this->col_indices != nullptr && this->vals != nullptr) && "Cannot scale: Matrix is empty.");
    std::transform(std::execution::par, this->vals->begin(), this->vals->end(), this->vals->begin(), [&factor](T val)
                   { return val * factor; });
}

// ============================= Helpers/Validators ==========================================================

template <typename T>
void CSRMatrix<T>::construct(const std::vector<T> &vals, const std::vector<size_t> &row_pointers, const std::vector<size_t> &col_indices) {
    this->row_pointers = new std::vector<size_t>(row_pointers);
    this->col_indices = new std::vector<size_t>(col_indices);
    this->vals = new std::vector<T>(vals);
}

template <typename T>
void CSRMatrix<T>::deepCopy(const CSRMatrix<T> &other) {
    assert((other.row_pointers != nullptr && other.col_indices != nullptr && other.vals != nullptr) && "Cannot copy: Input matrix is empty.");
    this->row_num = other.row_num;
    this->col_num = other.col_num;
    this->nnz = other.nnz;
    this->row_pointers = new std::vector<size_t>(*(other.row_pointers));
    this->col_indices = new std::vector<size_t>(*(other.col_indices));
    this->vals = new std::vector<T>(*(other.vals));
}

template <typename T>
void CSRMatrix<T>::validateCoordinates(size_t row, size_t col) const {
    // assert((row >= 0 && row < this->row_num) && "Row index out of bounds.");
    // assert((col >= 0 && col < this->col_num) && "Column index out of bounds.");

    if (row >= this->row_num)
        throw InvalidCoordinatesException("Row index out of bounds.");

    if (col >= this->col_num)
        throw InvalidCoordinatesException("Column index out of bounds.");
}

template <typename T>
void CSRMatrix<T>::validateCoordinates(size_t row) const {
    // assert((row >= 0 && row < this->row_num) && "Row index out of bounds.");

    if (row >= this->row_num)
        throw InvalidCoordinatesException("Row index out of bounds.");
}

template <typename T>
void CSRMatrix<T>::destruct() {
    delete this->vals;
    delete this->col_indices;
    delete this->row_pointers;
}

template <typename T>
void CSRMatrix<T>::shallowCopy(CSRMatrix<T> &&other) {
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
void CSRMatrix<T>::insert(size_t index, size_t row, size_t col, T val) {
    assert((this->row_pointers != nullptr && this->col_indices != nullptr && this->vals != nullptr) && "Cannot insert into an empty matrix.");
    this->vals->insert(this->vals->begin() + index, val);
    this->col_indices->insert(this->col_indices->begin() + index, col);
    this->nnz += 1;

    for (size_t i = row; i < this->row_num; ++i) {
        (*(this->row_pointers))[i + 1] += 1;
    }
}

template <typename T>
void CSRMatrix<T>::remove(size_t index, size_t row) {
    assert((this->row_pointers != nullptr && this->col_indices != nullptr && this->vals != nullptr) && "Cannot remove from an empty matrix.");
    this->vals->erase(this->vals->begin() + index);
    this->col_indices->erase(this->col_indices->begin() + index);
    this->nnz -= 1;

    for (size_t i = row; i < this->row_num; ++i) {
        (*(this->row_pointers))[i + 1] -= 1;
    }
}

template <typename T>
inline size_t CSRMatrix<T>::scanRowSize() {
    assert((this->row_pointers != nullptr) && "Row Pointers is null.");
    std::partial_sum(this->row_pointers->begin(), this->row_pointers->end(), this->row_pointers->begin());
    return (*(this->row_pointers))[this->row_num];
}

template <typename T>
inline bool CSRMatrix<T>::isEmpty() const {
    return (this->row_pointers == nullptr || this->col_indices == nullptr || this->vals == nullptr);
}

// ============================= Friend Functions ==========================================================
template <typename T>
bool operator==(const CSRMatrix<T> &lhs, const CSRMatrix<T> &rhs) {
    return (((lhs.vals == nullptr && rhs.vals == nullptr) ||
             (lhs.vals != nullptr && rhs.vals != nullptr && *(lhs.vals) == *(rhs.vals))) &&
            ((lhs.col_indices == nullptr && rhs.col_indices == nullptr) ||
             (lhs.col_indices != nullptr && rhs.col_indices != nullptr && *(lhs.col_indices) == *(rhs.col_indices))) &&
            (*(lhs.row_pointers) == *(rhs.row_pointers)));
}

template <typename T>
bool operator!=(const CSRMatrix<T> &lhs, const CSRMatrix<T> &rhs) {
    return !(lhs == rhs);
}

template <typename T>
std::ostream &operator<<(std::ostream &os, const CSRMatrix<T> &matrix) {
    assert((matrix.row_pointers != nullptr && matrix.col_indices != nullptr && matrix.vals != nullptr) && "Cannot print: Matrix is empty.");
    for (size_t i = 0; i < matrix.row_num; ++i) {
        for (size_t j = 0; j < matrix.col_num; ++j) {
            T val = matrix.get(i, j);

            if (val != T()) {
                os << "(" << i << ", " << j << "): ";
                os << val << "\n";
            }
        }
    }

    return os;
}

template <typename T>
void swap(CSRMatrix<T> &lhs, CSRMatrix<T> &rhs) {
    if (lhs != rhs)
    {
        std::swap(lhs.row_num, rhs.row_num);
        std::swap(lhs.col_num, rhs.col_num);
        std::swap(lhs.nnz, rhs.nnz);
        std::swap(lhs.row_pointers, rhs.row_pointers);
        std::swap(lhs.vals, rhs.vals);
        std::swap(lhs.col_indices, rhs.col_indices);
    }
}

template <typename T>
void spgemm_saad(const CSRMatrix<T> &A, const CSRMatrix<T> &B, CSRMatrix<T> &C);

template <typename T>
void spgemm_rmerge(const CSRMatrix<T> &A, const CSRMatrix<T> &B, CSRMatrix<T> &C);

template <typename Func>
void launchThreads(size_t row_num, Func&& f) {
    std::vector<std::thread> threads;
    size_t rows_per_thread = (row_num + num_threads - 1) / num_threads;

    for (int t = 0; t < num_threads; ++t) {
        size_t start = t * rows_per_thread;
        size_t end = std::min(row_num, start + rows_per_thread);
        threads.emplace_back(std::forward<Func>(f), start, end);
    }

    for (auto &thread : threads)
        thread.join();
}