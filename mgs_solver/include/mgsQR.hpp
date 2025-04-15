#pragma once

#include <algorithm>
#include <vector>
#include <cmath>
#include <numeric>

/* ----------------------- gramSchmidt ----------------------- */
/*  Given a matrix A of dimension m by n, this algorithm
    computes a QR decomposition of A, where Q is a unitary
    m by n matrix and R is a n by n upper triangular matrix
    and A = QR.

    Input variables:
        a   : 2d vector, the ith element of which should correspond to the
              ith column of the matrix A. During the algorithm, the columns
              of Q will replace the columns of A.
        r   : 2d vector  the ith column of the upper triangular matrix R
              will be stored in the ith index of r.
        m   : number of columns in A.
        n   : number of rows in A.
        full: TRUE  => full QR factorization computed
              FALSE => thin QR factorization computed

    Features: This implementation has time complexity O(m n^2)
    and requires O(1) additional memory.

    Remarks: Due to the nature of the problem, if A is nearly
    rank-deficient then the resulting columns of Q may not
    exhibit the orthogonality property.                        */

template <typename T>
void gramSchmidt(std::vector<std::vector<T>>& a, std::vector<T>& r, const size_t m, const size_t n)
{
    double anorm = 0.0;
}

template <typename T>
void mgsQRSolve(std::vector<std::vector<T>>& A, std::vector<T>& rhs, std::vector<T>& x, const size_t rows, const size_t cols)
{
    std::vector<std::vector<T>> R(cols, std::vector<T>(cols));
    gramSchmidt(A, R, rows, cols);
}