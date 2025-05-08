#include "CSRMatrixMock.hpp"
#include "testlib.hpp"

void _multiplicationFail1()
{
    // Create a matrix
    /*
            [0 2 0 1]
        A = [1 3 0 0]
            [0 7 4 0]
        row pointers:   [0, 2, 4, 6]
        column indices: [1, 3, 0, 1, 1, 2]
        values:         [2, 1, 1, 3, 7, 4]

    */
    std::vector<size_t> rowPointersA = {0, 2, 4, 6};
    std::vector<size_t> colIndicesA = {1, 3, 0, 1, 1, 2};
    std::vector<int> valsA = {2, 1, 1, 3, 7, 4};
    CSRMatrix<int> A(3, 4, valsA, rowPointersA, colIndicesA);
    std::vector<int> x(3, 1);
    A.multiply(x);
}

void testMultiplicationFail1()
{
    std::cout << "multiply() fail #1..." << std::flush;
    assertException("InvalidDimensionsException", _multiplicationFail1);
    std::cout << " OK" << std::endl;
}

void _multiplicationFail2()
{
    // Create two matrices with different dimensions
    /*
            [0 2 0 1]                                    [6 7 0 2 0 0]
        A = [1 3 0 0]                                    [0 1 0 0 0 0]
            [0 7 4 0]                               B =  [0 0 2 0 4 0]
        row pointers:   [0, 2, 4, 6]                     [7 3 2 0 0 0]
        column indices: [1, 3, 0, 1, 1, 2]               [1 1 0 0 3 4]
        values:         [2, 1, 1, 3, 7, 4]           row pointers:   [0, 3, 4, 6, 9, 13]
                                                     column indices: [0, 1, 3, 1, 2, 4, 0, 1, 2, 0, 1, 4, 5]
                                                     values:         [6, 7, 2, 1, 2, 4, 7, 3, 2, 1, 1, 3, 4]
    */
    std::vector<size_t> rowPointersA = {0, 2, 4, 6};
    std::vector<size_t> colIndicesA = {1, 3, 0, 1, 1, 2};
    std::vector<int> valsA = {2, 1, 1, 3, 7, 4};
    CSRMatrix<int> A(3, 4, valsA, rowPointersA, colIndicesA);

    std::vector<size_t> rowPointersB = {0, 3, 4, 6, 9, 13};
    std::vector<size_t> colIndicesB = {0, 1, 3, 1, 2, 4, 0, 1, 2, 0, 1, 4, 5};
    std::vector<int> valsB = {6, 7, 2, 1, 2, 4, 7, 3, 2, 1, 1, 3, 4};
    CSRMatrix<int> B(5, 6, valsB, rowPointersB, colIndicesB);

    A.multiply(B);
}

void testMultiplicationFail2()
{
    std::cout << "multiply() fail #2..." << std::flush;
    assertException("InvalidDimensionsException", _multiplicationFail2);
    std::cout << " OK" << std::endl;
}

void testVectorMultiplication()
{
    for (int N = 0; N < 5e3; ++N)
    {
        std::cout << "\rvector multiplication...#" << N + 1 << std::flush;

        // generate random matrices
        size_t rows = rand() % 16 + 1;
        size_t cols = rand() % 16 + 1;

        std::vector<int> vec = generateRandomVector<int>(cols);
        std::vector<std::vector<int>> classicMatrix = generateRandomMatrix<int>(rows, cols);
        CSRMatrixMock<int> sparseMatrix(classicMatrix);

        // calculate result manually
        std::vector<int> manualResult = multiplyMatrixByVector(classicMatrix, vec);

        // method
        assertEquals<std::vector<int>>(manualResult, *sparseMatrix.multiply(vec), "Incorrect matrix-vector multiplication result");

        // operator
        assertEquals<std::vector<int>>(manualResult, *(sparseMatrix * vec), "Incorrect matrix-vector multiplication result (operator *)");
    }

    std::cout << " OK" << std::endl;
}

void testMatrixMultiplication()
{
    for (int N = 0; N < 5e3; ++N)
    {
        std::cout << "\rmatrices multiplication... #" << N + 1 << std::flush;
        
        // generate random matrices
        size_t rowsA = rand() % 16 + 1;
        size_t colsArowsB = rand() % 16 + 1;
        size_t colsB = rand() % 16 + 1;

        std::vector<std::vector<int>> classicMatrixA = generateRandomMatrix<int>(rowsA, colsArowsB);
        CSRMatrixMock<int> sparseMatrixA(classicMatrixA);

        std::vector<std::vector<int>> classicMatrixB = generateRandomMatrix<int>(colsArowsB, colsB);
        CSRMatrixMock<int> sparseMatrixB(classicMatrixB);

        // calculate result manually
        std::vector<std::vector<int>> manualResult = multiplyMatrices<int>(classicMatrixA, classicMatrixB);

        // method
        assertEquals<CSRMatrix<int>, std::vector<std::vector<int>>>(
            *sparseMatrixA.multiply(sparseMatrixB),
            manualResult,
            "Incorrect matrix-matrix multiplication result");

        // operator
        assertEquals<CSRMatrix<int>, std::vector<std::vector<int>>>(
            *(sparseMatrixA * sparseMatrixB),
            manualResult,
            "Incorrect matrix-matrix multiplication result (operator *)");
    }

    std::cout << " OK" << std::endl;
}