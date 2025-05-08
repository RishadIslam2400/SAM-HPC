#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "CSRMatrix.hpp"

/**
 * @brief: Read a sparse matrix stored in CSR format in a text file. The first line in the
 * text file contains the number of rows, columns, and non-zero elements. The following lines contain
 * the row pointers, column indices, and values of the non-zero elements. The function populates the
 * the provided CSRMatrix object with the data read from the file.
 *
 * @param filename: A C-style string representing the path to the input text file containing the matrix.
 * @param matrix: A pointer to an uninitialized or empty CSRMatrix<double> object that will be filled
 *                with the matrix data read from the file.
 * @return true If the matrix was read successfully, false otherwise.
 */
template <typename T>
bool read_mat(const char *filename, CSRMatrix<T> *matrix)
{
    std::ifstream infile(filename);
    if (!infile.is_open())
    {
        std::cerr << "Error opening file: " << filename << std::endl;
        return false;
    }

    size_t rows, cols, nnz;
    std::string line;

    // First line: rows, cols, nnz
    if (!std::getline(infile, line))
    {
        std::cerr << "Error reading matrix dimensions from file: " << filename << std::endl;
        return false;
    }

    std::istringstream header(line);
    if (!(header >> rows >> cols >> nnz))
    {
        std::cerr << "Error parsing matrix dimensions from file: " << filename << std::endl;
        return false;
    }

    matrix->row_num = rows;
    matrix->col_num = cols;
    matrix->nnz = nnz;

    // Second line: row pointers
    if (!std::getline(infile, line))
    {
        std::cerr << "Error reading row pointers from file: " << filename << std::endl;
        return false;
    }
    std::istringstream rowPtrStream(line);
    matrix->row_pointers = new std::vector<size_t>(rows + 1);
    for (size_t i = 0; i <= rows; ++i)
    {
        if (!(rowPtrStream >> (*(matrix->row_pointers))[i]))
        {
            std::cerr << "Error parsing row pointers from file: " << filename << std::endl;
            return false;
        }
    }

    // Third line: column indices
    if (!std::getline(infile, line))
    {
        std::cerr << "Error reading column indices from file: " << filename << std::endl;
        return false;
    }
    std::istringstream colIdxStream(line);
    matrix->col_indices = new std::vector<size_t>(nnz);
    for (size_t i = 0; i < nnz; ++i)
    {
        if (!(colIdxStream >> (*(matrix->col_indices))[i]))
        {
            std::cerr << "Error parsing column indices from file: " << filename << std::endl;
            return false;
        }
    }

    // Fourth line: values
    if (!std::getline(infile, line))
    {
        std::cerr << "Error reading values from file: " << filename << std::endl;
        return false;
    }
    std::istringstream valStream(line);
    matrix->vals = new std::vector<T>(nnz);
    for (size_t i = 0; i < nnz; ++i)
    {
        if (!(valStream >> (*(matrix->vals))[i]))
        {
            std::cerr << "Error parsing values from file: " << filename << std::endl;
            return false;
        }
    }

    infile.close();
    return true;
}

/**
 * @brief: Read a sparse matrix stored in CSR format in a text file. The first line in the
 * text file contains the number of rows, columns, and non-zero elements. The following lines contain
 * the row pointers, column indices, and values of the non-zero elements. The function creates a new CSRMatrix object 
 * with the data read from the file.
 *
 * @param filename: A C-style string representing the path to the input text file containing the matrix.
 * @param matrix: A pointer to an uninitialized or empty CSRMatrix<double> object that will be filled
 *                with the matrix data read from the file.
 * @return A CSRMatrix object containing the matrix data read from the file.
 */
template <typename T>
CSRMatrix<T>  read_mat(const char *filename)
{
    std::ifstream infile(filename);
    if (!infile.is_open())
    {
        std::cerr << "Error opening file: " << filename << std::endl;
        return CSRMatrix<T>();
    }

    size_t rows, cols, nnz;
    std::string line;

    // First line: rows, cols, nnz
    if (!std::getline(infile, line))
    {
        std::cerr << "Error reading matrix dimensions from file: " << filename << std::endl;
        return CSRMatrix<T>();
    }

    std::istringstream header(line);
    if (!(header >> rows >> cols >> nnz))
    {
        std::cerr << "Error parsing matrix dimensions from file: " << filename << std::endl;
        return CSRMatrix<T>();
    }

    // Second line: row pointers
    if (!std::getline(infile, line))
    {
        std::cerr << "Error reading row pointers from file: " << filename << std::endl;
        return CSRMatrix<T>();
    }
    std::istringstream rowPtrStream(line);
    std::vector<size_t> row_pointers(rows + 1);
    for (size_t i = 0; i <= rows; ++i)
    {
        if (!(rowPtrStream >> row_pointers[i]))
        {
            std::cerr << "Error parsing row pointers from file: " << filename << std::endl;
            return CSRMatrix<T>();
        }
    }

    // Third line: column indices
    if (!std::getline(infile, line))
    {
        std::cerr << "Error reading column indices from file: " << filename << std::endl;
        return CSRMatrix<T>();
    }
    std::istringstream colIdxStream(line);
    std::vector<size_t> col_indices(nnz);
    for (size_t i = 0; i < nnz; ++i)
    {
        if (!(colIdxStream >> col_indices[i]))
        {
            std::cerr << "Error parsing column indices from file: " << filename << std::endl;
            return CSRMatrix<T>();
        }
    }

    // Fourth line: values
    if (!std::getline(infile, line))
    {
        std::cerr << "Error reading values from file: " << filename << std::endl;
        return CSRMatrix<T>();
    }
    std::istringstream valStream(line);
    std::vector<T> vals(nnz);
    for (size_t i = 0; i < nnz; ++i)
    {
        if (!(valStream >> vals[i]))
        {
            std::cerr << "Error parsing values from file: " << filename << std::endl;
            return CSRMatrix<T>();
        }
    }

    infile.close();

    return CSRMatrix<T>(rows, cols, nnz, vals, row_pointers, col_indices);
}