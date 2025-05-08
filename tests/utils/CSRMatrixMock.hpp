#pragma once

#include <vector>

#include "CSRMatrix.hpp"
#include "spgemm.hpp"

/**
 * This class is used only for testing purposes
 *
 * @internal
 */
template <typename T>
class CSRMatrixMock : public CSRMatrix<T>
{
public:
    CSRMatrixMock(const size_t rows, const size_t cols, const std::vector<T> &vals,
                  const std::vector<size_t> &row_pointers, const std::vector<size_t> &col_indices)
        : CSRMatrix<T>(rows, cols, vals, row_pointers, col_indices) {}
    CSRMatrixMock(const size_t rows, const size_t cols, const size_t nnz, const std::vector<T> &vals,
                  const std::vector<size_t> &row_pointers, const std::vector<size_t> &col_indices)
        : CSRMatrix<T>(rows, cols, nnz, vals, row_pointers, col_indices) {}
    CSRMatrixMock(const std::vector<T> &vals, const std::vector<size_t> &row_pointers, const std::vector<size_t> &col_indices)
        : CSRMatrix<T>(vals, row_pointers, col_indices) {}
    CSRMatrixMock(const size_t rows, const size_t cols, const T *vals,
                  const size_t *row_pointers, const size_t *col_indices)
        : CSRMatrix<T>(rows, cols, vals, row_pointers, col_indices) {}
    CSRMatrixMock(const size_t rows, const size_t cols, const size_t nnz, const T *vals,
                  const size_t *row_pointers, const size_t *col_indices)
        : CSRMatrix<T>(rows, cols, nnz, vals, row_pointers, col_indices) {}
    CSRMatrixMock(const std::vector<std::vector<T>> &matrix) : CSRMatrix<T>(matrix) {}

    /**
     * Sends internal storage info to given output stream
     *
     * @param os output stream
     * @return void
     */
    void printInfo(std::ostream &os) const
    {
        os << "rows (" << (this->row_pointers->size() - 1) << "): [";
        for (std::vector<size_t>::const_iterator intIt = this->row_pointers->cbegin(); intIt != this->row_pointers->cend(); ++intIt)
        {
            if (intIt > this->row_pointers->cbegin())
            {
                os << ", ";
            }

            os << *intIt;
        }
        os << "]" << std::endl;

        os << "cols";
        if (this->col_indices == nullptr)
        {
            os << ": NULL" << std::endl;
        }
        else
        {
            os << " (" << this->col_indices->size() << "): [";
            for (std::vector<size_t>::const_iterator intIt = this->col_indices->cbegin(); intIt != this->col_indices->cend(); ++intIt)
            {
                if (intIt > this->col_indices->cbegin())
                {
                    os << ", ";
                }

                os << *intIt;
            }
            os << "]" << std::endl;
        }

        os << "vals";
        if (this->vals == nullptr)
        {
            os << ": NULL" << std::endl;
        }
        else
        {
            os << " (" << this->vals->size() << "): [";
            for (typename std::vector<T>::const_iterator valIt = this->vals->cbegin(); valIt != this->vals->cend(); ++valIt)
            {
                if (valIt > this->vals->cbegin())
                {
                    os << ", ";
                }

                os << *valIt;
            }
            os << "]" << std::endl;
        }
    }
};

template <typename T>
bool operator==(const CSRMatrix<T> &sparse, const std::vector<std::vector<T>> &classical)
{
    for (size_t i = 0; i < classical.size(); ++i)
    {
        for (size_t j = 0; j < classical[i].size(); ++j)
        {
            if (sparse.get(i, j) != classical[i][j])
                return false;
        }
    }

    return true;
}