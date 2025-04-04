#pragma once

#include <vector>
#include <iostream>
#include <concepts>

template <typename T>
    requires std::is_arithmetic_v<T>
std::vector<T> generateRandomVector(int size)
{
    std::vector<T> vector(size, T());

    for (int i = 0; i < size; ++i)
    {
        vector[i] = rand() % 101;
    }

    return vector;
}

template <typename T>
    requires std::is_arithmetic_v<T>
std::vector<std::vector<T>> generateRandomMatrix(int rows, int cols)
{
    std::vector<std::vector<T>> matrix(rows, std::vector<T>(cols, T()));

    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
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
    int rows = a.size();
    int cols = a.front().size();
    std::vector<std::vector<T>> result(rows, std::vector<T>(cols, T()));

    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
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
    int rows = a.size();
    int cols = a.front().size();
    std::vector<std::vector<T>> result(rows, std::vector<T>(cols, T()));

    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
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
    int rows = static_cast<int>(matrix.size());
    int cols = static_cast<int>(matrix[0].size());

    std::vector<T> result(rows, T());

    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
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
    int rowsA = a.size();
    int colsA = a.front().size();
    int colsB = b.front().size();

    std::vector<std::vector<T>> result(rowsA, std::vector<T>(colsB, T()));

    for (int i = 0; i < rowsA; ++i)
    {
        for (int j = 0; j < colsB; ++j)
        {
            result[i][j] = T();
            for (int k = 0; k < colsA; ++k)
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
    for (int i = 0; i < static_cast<int>(vec.size()); ++i)
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