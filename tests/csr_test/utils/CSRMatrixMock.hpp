#pragma once

#include <vector>

#include "CSRMatrix.hpp"

/**
 * This class is used only for testing purposes
 *
 * @internal
 */
template <typename T>
class CSRMatrixMock : public SparseMatrix::CSRMatrix<T>
{
public:
    CSRMatrixMock(int n) : SparseMatrix::CSRMatrix<T>(n) {}
    CSRMatrixMock(int rows, int cols) : SparseMatrix::CSRMatrix<T>(rows, cols) {}
    CSRMatrixMock(int rows, int cols, const std::vector<T> &vals, const std::vector<int> &row_pointers, const std::vector<int> &col_indices)
        : SparseMatrix::CSRMatrix<T>(rows, cols, vals, row_pointers, col_indices) {}

    /**
     * Sends internal storage info to given output stream
     *
     * @param os output stream
     * @return void
     */
    void printInfo(std::ostream &os) const
    {
        os << "rows (" << (this->row_pointers->size() - 1) << "): [";
        for (std::vector<int>::const_iterator intIt = this->row_pointers->cbegin(); intIt != this->row_pointers->cend(); ++intIt)
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
            for (std::vector<int>::const_iterator intIt = this->col_indices->cbegin(); intIt != this->col_indices->cend(); ++intIt)
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

    static const CSRMatrixMock<T> fromVectors(const std::vector<std::vector<T>>& vec)
    {
        const int row_size = static_cast<int>(vec.size());
        const int col_size = static_cast<int>(vec[0].size());
        CSRMatrixMock<T> matrix(row_size, col_size);

        for (int i = 0; i < row_size; ++i)
        {
            for (int j = 0; j < col_size; ++j)
            {
                matrix.set(vec[i][j], i, j);
            }
        }

        return matrix;
    }
};

template <typename T>
bool operator==(const SparseMatrix::CSRMatrix<T> &sparse, const std::vector<std::vector<T>> &classical)
{
    for (int i = 0; i < static_cast<int>(classical.size()); ++i)
    {
        for (int j = 0; j < static_cast<int>(classical[i].size()); ++j)
        {
            if (sparse.get(i, j) != classical[i][j])
                return false;
        }
    }

    return true;
}