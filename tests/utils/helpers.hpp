#pragma once

#include <vector>
#include <iostream>
#include <concepts>

template <typename T>
    requires std::is_arithmetic_v<T>
std::vector<T> generateRandomVector(size_t size)
{
    std::vector<T> vector(size, T());

    for (size_t i = 0; i < size; ++i)
    {
        vector[i] = rand() % 101;
    }

    return vector;
}

template <typename T>
    requires std::is_arithmetic_v<T>
std::vector<std::vector<T>> generateRandomMatrix(size_t rows, size_t cols)
{
    std::vector<std::vector<T>> matrix(rows, std::vector<T>(cols, T()));

    for (size_t i = 0; i < rows; ++i)
    {
        for (size_t j = 0; j < cols; ++j)
        {
            matrix[i][j] = rand() % 101;
        }
    }

    return matrix;
}

template <typename T>
    requires std::is_arithmetic_v<T>
std::vector<std::vector<T>> addMatrices(const std::vector<std::vector<T>> &a, const std::vector<std::vector<T>> &b)
{
    size_t rows = a.size();
    size_t cols = a.front().size();
    std::vector<std::vector<T>> result(rows, std::vector<T>(cols, T()));

    for (size_t i = 0; i < rows; ++i)
    {
        for (size_t j = 0; j < cols; ++j)
        {
            result[i][j] = a[i][j] + b[i][j];
        }
    }

    return result;
}

template <typename T>
    requires std::is_arithmetic_v<T>
std::vector<std::vector<T>> subtractMatrices(const std::vector<std::vector<T>> &a, const std::vector<std::vector<T>> &b)
{
    size_t rows = a.size();
    size_t cols = a.front().size();
    std::vector<std::vector<T>> result(rows, std::vector<T>(cols, T()));

    for (size_t i = 0; i < rows; ++i)
    {
        for (size_t j = 0; j < cols; ++j)
        {
            result[i][j] = a[i][j] - b[i][j];
        }
    }

    return result;
}

template <typename T>
    requires std::is_arithmetic_v<T>
std::vector<T> multiplyMatrixByVector(const std::vector<std::vector<T>> &matrix, const std::vector<T> &vec)
{
    size_t rows = matrix.size();
    size_t cols = matrix[0].size();

    std::vector<T> result(rows, T());

    for (size_t i = 0; i < rows; ++i)
    {
        for (size_t j = 0; j < cols; ++j)
        {
            result[i] += matrix[i][j] * vec[j];
        }
    }
    return result;
}

template <typename T>
    requires std::is_arithmetic_v<T>
std::vector<std::vector<T>> multiplyMatrices(const std::vector<std::vector<T>> &a, const std::vector<std::vector<T>> &b)
{
    size_t rowsA = a.size();
    size_t colsA = a.front().size();
    size_t colsB = b.front().size();

    std::vector<std::vector<T>> result(rowsA, std::vector<T>(colsB, T()));

    for (size_t i = 0; i < rowsA; ++i)
    {
        for (size_t j = 0; j < colsB; ++j)
        {
            result[i][j] = T();
            for (size_t k = 0; k < colsA; ++k)
            {
                result[i][j] += a[i][k] * b[k][j];
            }
        }
    }

    return result;
}

template <typename T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &vec)
{
    os << "[";
    for (size_t i = 0; i < vec.size(); ++i)
    {
        if (i > 0)
        {
            os << ", ";
        }
        os << vec[i];
    }
    os << "]";
    return os;
}

bool operator==(const std::vector<double> &lhs, const std::vector<double> &rhs)
{
    if (lhs.size() != rhs.size())
        return false;
    for (size_t i = 0; i < lhs.size(); ++i)
    {
        if (std::abs(lhs[i] - rhs[i]) > 1e-7)
            return false;
    }
    return true;
}