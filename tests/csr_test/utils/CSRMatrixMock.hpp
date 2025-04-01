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

};