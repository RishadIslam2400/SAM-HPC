#pragma once

#include <algorithm>
#include <vector>
#include <cmath>
#include <numeric>
#include <cassert>
#include <span>

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
void gramSchmidt(std::vector<std::vector<double>> &a, std::vector<std::vector<double>> &q, std::vector<std::vector<double>> &r, const size_t m, const size_t n)
{
    double anorm = 0.0;
    const double tol = 10e-10;
    for (size_t i = 0; i < n; ++i)
    {
        // Compute the norm of the ith column of A
        // r_ii = ||a_i||
        r[i][i] = std::sqrt(std::transform_reduce(a[i].begin(), a[i].end(), 0.0, std::plus<double>(), [](double x)
                                                  { return x * x; }));

        // Normalize the ith column of A
        // a_i = a_i / r_ii
        if (r[i][i] > tol)
        {
            std::transform(a[i].begin(), a[i].end(), a[i].begin(), [&](double x)
                           { return x / r[i][i]; });
        }
        else if (i == 0)
        {
            // set the first column a[0] = [1, 0, 0, ...]^T
            a[i][0] = 1.0;
            for (size_t j = 1; j < m; ++j)
            {
                a[i][j] = 0.0;
            }
        }
        else
        {
            // Create a_i orthogonal to the previous columns
            double scale_factor = -a[0][i];
            std::transform(a[i].begin(), a[i].end(), a[i].begin(), [&](double x)
                           { return scale_factor * x; });
            a[i][i] += 1.0;
            for (size_t j = 1; j < i; ++j)
            {
                std::transform(a[j].begin(), a[j].end(), a[i].begin(), a[i].begin(), [&](double x, double y)
                               { return y - x * a[j][i]; });
            }
            anorm = std::sqrt(std::transform_reduce(a[i].begin(), a[i].end(), 0.0, std::plus<double>(), [](double x)
                                                    { return x * x; }));            
            std::transform(a[i].begin(), a[i].end(), a[i].begin(), [&](double x)
                           { return x / anorm; });
        }

        for (size_t j = i + 1; j < n; ++j)
        {
            // Compute the inner product of the ith and jth columns
            // r_ij = a_i * a_j
            r[j][i] = std::transform_reduce(a[i].begin(), a[i].end(), a[j].begin(), 0.0, std::plus<double>(), [&](double x, double y)
                                            { return x * y; });

            // a_j -= r_ij * a_i
            std::transform(a[i].begin(), a[i].end(), a[j].begin(), a[j].begin(), [&](double xi, double yi)
                           { return yi - r[j][i] * xi; });
        }
    }
}


void mgsQRSolve(std::vector<std::vector<double>> &A, std::vector<double> &rhs, std::span<double> x, const size_t cols, const size_t rows)
{
    if (rows == 0 || cols == 0)
    {
        return;
    }
    
    assert((A.size() == cols && A[0].size() == rows) && "Matrix A must be of size rows x cols.");
    assert(x.size() == cols && "Vector X must be equal to the number of columns.");
    assert(rhs.size() == rows && "RHS vector must be equal to the number of rows.");
    // @todo: handle underdetermined systems
    assert(rows >= cols && "Can't solve underdetermined linear system.");
    
    std::vector<std::vector<double>> R(cols, std::vector<double>(cols));
    std::vector<std::vector<double>> Q(rows, std::vector<double>(rows));

    gramSchmidt(A, Q, R, rows, cols);

    // Solve the system R * x = Q^T * b
    std::vector<double> QTb(cols, 0.0);

    // Compute Q^T * b
    // Accesing the rows in a column-major matrix
    for (size_t i = 0; i < rows; ++i)
    {
        for (size_t j = 0; j < cols; ++j)
        {
            QTb[j] += A[j][i] * rhs[i];
        }
    }

    // Back substitution to solve R * x = Q^T * b
    for (int i = static_cast<int>(cols) - 1; i >= 0; --i)
    {
        x[i] = QTb[i];
        for (size_t j = i + 1; j < cols; ++j)
        {
            x[i] -= R[j][i] * x[j];
        }

        x[i] /= R[i][i]; // Avoid division by zero
    }
}