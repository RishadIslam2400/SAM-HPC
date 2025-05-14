#pragma once

#include <algorithm>
#include <vector>
#include <cmath>
#include <numeric>
#include <cassert>
#include <span>

/* ------------------ partialdot_product ------------------ */
/*  Given two vectors of the same length and an index this function returns
    the value of the dot product of the two vectors

    Input variables:
        x     : first vector.
        y     : second vector.
        index : starting index for the vector

    Example: Suppose x is a vector {1, 2, 3, 4}, y is a vector {5, 6, 7, 8} and index = 1
    Then the value returned by partialdot_product(x, y) is 65, which is computed by
        x[1] * y[1] + x[2] * y[2] + x[3] * y[3]
        = 2 * 6 + 3 * 7 + 4 * 8
        = 12 + 21 + 32
        = 65

    Features: This implementation has time complexity
    O(x.size()) and requires O(1) additional memory.          */
double partial_dot_product(const std::vector<double> &x, const std::vector<double> &y, const size_t index) {
    assert(x.size() == y.size() && "Vectors must be of the same size.");
    assert(index <= x.size() && "Index out of bounds.");

    return std::transform_reduce(x.begin() + index, x.end(), y.begin() + index, 0.0, std::plus<double>(), std::multiplies<double>());
}

/* ----------------------- scalar_div ----------------------- */
/*  Given a vector, and a scalar value this function divides the values from the
    vector by the scalar value and stores the computed number in the same vector.

    Input variables:
        x     : vector whose components are to be divided by r
        r     : scalar used in division.


    Features: This implementation has time complexity
    O(length) and requires O(1) additional memory.            */

void scalar_div(std::vector<double> &x, const double r) {
    std::transform(x.begin(), x.end(), x.begin(), [&r](double a) { return a / (r + 1e-13); });
}

/* -------------------- subdot_product -------------------- */
/*  Given two vectors and an index this function returns
    the dot product of the two vectors x[index : index + length] and
    y[0 : length]. It is necessary that index + length is
    at most the length of the first array.

    Input variables:
        x     : first input vector.
        y     : second input vector.
        index : starting index for subarray of x.

    Example: Suppose x is a vector {1, 2, 3, 4, 5}, y is a vector {-1, -2, -3},
    and index = 2. Then the value returned by executing subdot_product(x, y, 2)
    is -26, which is computed by
    x[2] * y[0] + x[3] * y[1] + x[4] * y[2]

    =  3   *  -1  +  4   *  -2  +  5   *  -3

    = -    3      -      8      -      15

    = -26.

    Features: This implementation has time complexity
    O(length) and requires O(1) additional memory.          */
double subdot_product(const std::vector<double> &x, const std::vector<double> &y, const size_t index) {
    assert(index <= x.size() && "Index out of bounds.");
    assert(index + y.size() <= x.size() && "Vectors must be of the same size.");

    return std::transform_reduce(x.begin() + index, x.end(), y.begin(), 0.0, std::plus<double>(), [](double a, double b) { return a * b; });
}

/* --------------------- partialscalar_sub --------------------- */
/*  Given two vectors, a scalar value, and an index, this function
    multiplies the the first vector by the scalar value and then
    subtracts the computed components from the components the second
    vector starting at the given index.

    Input variables:
        x     : vector whose components are to be multiplied by r then
                subtracted from the components of the second array, y.
        r     : scalar used in multiplication.
        index :
        y     : vector in which the components of x are to be stored.

    Example: Suppose x is a vector {3, 4, 5}, y is a vector {1, 2, 0, 0, 0},
    r = -1, and index = 2. Then after executing partialscalar_sub(x, -1, 2, y),
    the array pointed to by y is now {1, 2, -3, -4, -5}.

    Features: This implementation has time complexity
    O(length) and requires O(1) additional memory.               */
void partialscalar_sub(const std::vector<double> &x, const double r, const size_t index, std::vector<double> &y) {
    assert(index <= y.size() && "Index out of bounds.");
    assert(index + x.size() <= y.size() && "Vectors must be of the same size.");

    std::transform(x.begin(), x.end(), y.begin() + index, y.begin() + index, [r](double a, double b) { return b - r * a; });
}

/* ----------------------- householder ----------------------- */
/*  Given a matrix A of dimension m by n (with n <= m) and
    arrays v_i of dimension m-i, for i = 1, ..., m - 1,
    respectively, this algorithm computes n reflection vectors
    and the factor R of a full QR decomposition of A, where R
    is a m by n upper triangular matrix. The n reflection
    vectors are stored in the arrays v_1, ..., v_n and the
    columns of A are overwritten by the columns of R.

    Input variables:
        a: reference to a 2d vector, the ith element of
            which should correspond to the ith column of the
            matrix A. During the algorithm, the columns of R
            will overwrite the columns of A.
        v: reference to vector of vectors in which the ith
            reflection vector of dimension m - i will be
            stored.
        m: number of rows of A.
        n: number of columns of A.

    Features: The number of flops for this implementation is
    ~ 2 * m * n^2 - (2/3) * n^3 and requires O(1) additional
    memory.                                                    */

void householderQR(std::vector<std::vector<double>> &a, std::vector<std::vector<double>> &v, const size_t m, const size_t n) {
    for (size_t i = 0; i < n; ++i) {
        // Set v[i] equal to subvector a[i][i:m]
        v[i].assign(a[i].begin() + i, a[i].end());
        

        /* vpartdot = ||v[i]||^2 - v[i][0] * v[i][0]; since vpartdot
           is unaffected by the change in v[i][0], storing this value
           prevents the need to recalculate the entire norm of v[i]
           after updating v[i][0] in the following step              */
        const double vpartdot = partial_dot_product(v[i], v[i], 1);

        // set v[i][0] = v[i][0] + sign(v[i][0]) * ||v[i]||;
        if (v[i][0] < 0) {
            v[i][0] -= std::sqrt(v[i][0] * v[i][0] + vpartdot);
        } else {
            v[i][0] += std::sqrt(v[i][0] * v[i][0] + vpartdot);
        }

        // Normalize v[i]
        const double vnorm = std::sqrt(v[i][0] * v[i][0] + vpartdot);
        scalar_div(v[i], vnorm);

        for (size_t j = i; j < n; ++j) {
            // Set a[j][i:m] = a[j][i:m] - 2 * v[i] * (v[i]^T * a[j][i:m])
            double vTa = subdot_product(a[j], v[i], i);
            vTa *= 2.0;
            partialscalar_sub(v[i], vTa, i, a[j]);
        }
    }
}

void householderQRSolve(std::vector<std::vector<double>> &A, std::vector<double> &rhs, std::span<double> x, const size_t cols, const size_t rows) {
    if (rows == 0 || cols == 0) {
        return;
    }

    assert((A.size() == cols && A[0].size() == rows) && "Matrix A must be of size rows x cols.");
    assert(x.size() == cols && "Vector X must be equal to the number of columns.");
    assert(rhs.size() == rows && "RHS vector must be equal to the number of rows.");
    // @todo: handle underdetermined systems
    assert(rows >= cols && "Can't solve underdetermined linear system.");

    // Initialize the reflector vectors
    std::vector<std::vector<double>> V(cols);
    householderQR(A, V, rows, cols);

    // Apply Q^T to the right-hand side
    for (size_t i = 0; i < cols; ++i) {
        // Compute 2 * (V[i]^T * rhs[i:rows]) * V[i]
        const double vTy = subdot_product(rhs, V[i], i) * 2.0;
        partialscalar_sub(V[i], vTy, i, rhs);
    }

    // Solve Rx = Q^T b using back-substitution
    for (int i = static_cast<int>(cols) - 1; i >= 0; --i) {
        x[i] = rhs[i];
        for (size_t j = i + 1; j < cols; ++j)
        {
            x[i] -= A[j][i] * x[j];
        }
        x[i] /= A[i][i];
    }
}