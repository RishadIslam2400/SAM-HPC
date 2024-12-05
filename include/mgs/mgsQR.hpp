#pragma once

#include <algorithm>
#include <vector>
#include <cmath>
#include <execution>
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

void gramSchmidt(std::vector<std::vector<double>>& a, std::vector<std::vector<double>>& r, size_t m, size_t n) {
    double anorm{0.0};
    double tol = 10e-7; // Paramter can be passed

    for (size_t i = 0; i < n; ++i) {
        // Compute the norm of the ith column of A
        // r_ii = ||a_i||
        r[i][i] = std::sqrt(std::transform_reduce(
            std::execution::seq,
            a[i].begin(),
            a[i].end(),
            0.0,
            std::plus<double>(),
            [](double val) { return val * val; }
        ));

        /* double sum = 0.0;
        for (size_t j = 0; j < m; ++j) {
            sum += a[i][j] * a[i][j];
        }
        r[i][i] = std::sqrt(sum); */

        // Normalize the ith column of A
        // a_i = a_i / ||a_i||
        if (r[i][i] > tol) {
            std::transform(
                std::execution::seq,
                a[i].begin(),
                a[i].end(),
                a[i].begin(),
                [&](double& elem) { return elem / r[i][i]; }
            );

            /* for (size_t j = 0; j < m; ++j) {
                a[i][j] /= r[i][i];
            } */
        } else if (i == 0) { // set a[0] = [1 0 0 ... 0]^T
            a[i][0] = 1;
            for (size_t j = 1; j < m; ++j) {
                a[i][j] = 0;
            }
        } else { // need to choose a_i othogonal to <a_1, a_2, ... a_i-1>
            double scale_factor = -a[0][i];
            std::transform(
                std::execution::seq,
                a[i].begin(),
                a[i].end(),
                a[i].begin(),
                [&](double& elem) { return scale_factor * elem; }
            );

            /* for (size_t j = 0; j < m; ++j) {
                a[i][j] = scale_factor * a[0][j];
            } */

            a[i][i] += 1;

            for (size_t j = 1; j < i; ++j) {
                std::transform(
                    std::execution::seq,
                    a[j].begin(),
                    a[j].end(),
                    a[i].begin(),
                    a[i].begin(),
                    [&](double xi, double yi) { return yi - xi * a[j][i]; }
                );

                /* for (size_t k = 0; k < m; ++k) {
                    a[i][k] -= a[j][k] * a[j][i];
                } */
            }

            anorm = std::sqrt(std::transform_reduce(
                                std::execution::seq,
                                a[i].begin(),
                                a[i].end(),
                                0.0,
                                std::plus<double>(),
                                [](double val) { return val * val; }
                            ));

            /*
            for (size_t j = 0; j < m; ++j) {
                anorm += a[i][j] * a[i][j];
            }
            anorm = std::sqrt(anorm); */
            
            std::transform(
                std::execution::seq,
                a[i].begin(),
                a[i].end(),
                a[i].begin(),
                [&](double& elem) { return elem / anorm; }
            );

            /* for (size_t j = 0; j < m; ++j) {
                a[i][j] /= anorm;
            } */
        }

        for (size_t j = i + 1; j < n; ++j) {
            // r_ij = a_i*a_j
            r[j][i] = std::transform_reduce(
                std::execution::seq,
                a[i].begin(),
                a[i].end(),
                a[j].begin(),
                0.0,
                std::plus<double>(),
                [](double xi, double yi) { return xi * yi; }
            );

            /* for (size_t k = 0; k < m; ++k) {
                r[j][i] += a[i][k] * a[j][k];
            } */

            // a_j -= r_ij * a_i
            std::transform(
                    std::execution::seq,
                    a[i].begin(),
                    a[i].end(),
                    a[j].begin(),
                    a[j].begin(),
                    [&](double xi, double yi) { return yi - r[j][i] * xi; }
            );

            /* for (size_t k = 0; k < m; ++k) {
                a[j][k] -= r[j][i] * a[i][k];
            } */
        }
    }
}

void mgsQRSolve(std::vector<std::vector<double>>& A, std::vector<double>& rhs, std::vector<double>& x, size_t rowDim, size_t colDim) {
    std::vector<std::vector<double>> R(colDim, std::vector<double>(colDim, 0.0));

    // Perfrom the QR decomposition
    gramSchmidt(A, R, rowDim, colDim);

    // Compute Q^T * b
    std::vector<double> QTb(colDim, 0.0);
    // Accessing the row in a column major matrix is not effecient (for large matrices)
    // Alternative: Tranpose the matrix before computation. Adds extra computational overhead.
    for (size_t i = 0; i < rowDim; ++i) {
        for (size_t j = 0; j < colDim; ++j) {
            QTb[j] += A[j][i] * rhs[i]; // Dot product of i-th row of Q with b
        }
    }

    // Back substitution to solve R * x = Q^T * b
    for (int i = static_cast<int>(colDim) - 1; i >= 0; --i) {
        x[i] = QTb[i];
        for (size_t j = i + 1; j < colDim; ++j) {
            x[i] -= R[j][i] * x[j];
        }

        x[i] /= R[i][i];
    }
}