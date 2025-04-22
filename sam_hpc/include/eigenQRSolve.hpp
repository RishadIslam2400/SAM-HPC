#pragma once

#include <eigen3/Eigen/Dense>

void eigenQRSolve(std::vector<std::vector<double>> &A, std::vector<double> &rhs, std::vector<double> &x, const size_t rows, const size_t cols)
{
    assert((A.size() == cols && A[0].size() == rows) && "Matrix A must be of size cols x rows.");
    assert(x.size() == cols && "Vector X must be equal to the number of columns.");
    assert(rhs.size() == rows && "RHS vector must be equal to the number of rows.");

    // Convert A and rhs to Eigen matrices
    Eigen::MatrixXd eigenA(rows, cols);
    Eigen::VectorXd eigenRhs(rows);
    for (size_t i = 0; i < cols; ++i)
    {
        for (size_t j = 0; j < rows; ++j)
        {
            eigenA(i, j) = A[i][j];
        }
    }
    for (size_t i = 0; i < rows; ++i)
    {
        eigenRhs(i) = rhs[i];
    }

    // Compute the QR decomposition of A
    Eigen::VectorXd eigenSol = eigenA.colPivHouseholderQr().solve(eigenRhs);

    // Convert the solution back to std::vector
    for (size_t i = 0; i < cols; ++i)
    {
        x[i] = eigenSol(i);
    }
}