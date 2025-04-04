#include "CSRMatrixMock.hpp"
#include "testlib.hpp"

void _multiplicationFail1()
{
    SparseMatrix::CSRMatrix<int> m(3, 4);
    std::vector<int> x(3, 1);
    m.multiply(x);
}

void testMultiplicationFail1()
{
    std::cout << "multiply() fail #1..." << std::flush;
    assertException("InvalidDimensionsException", _multiplicationFail1);
    std::cout << " OK" << std::endl;
}

void _multiplicationFail2()
{
    SparseMatrix::CSRMatrix<int> a(3, 4), b(5, 6);
    a.multiply(b);
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
        int rows = rand() % 16 + 1;
        int cols = rand() % 16 + 1;

        std::vector<int> vec = generateRandomVector<int>(cols);
        std::vector<std::vector<int>> classicMatrix = generateRandomMatrix<int>(rows, cols);
        CSRMatrixMock<int> sparseMatrix = CSRMatrixMock<int>::fromVectors(classicMatrix);

        // calculate result manually
        std::vector<int> manualResult = multiplyMatrixByVector(classicMatrix, vec);

        // method
        assertEquals<std::vector<int>>(manualResult, sparseMatrix.multiply(vec), "Incorrect matrix-vector multiplication result");

        // operator
        assertEquals<std::vector<int>>(manualResult, sparseMatrix * vec, "Incorrect matrix-vector multiplication result (operator *)");
    }

    std::cout << " OK" << std::endl;
}

void testMatrixMultiplication()
{
    for (int N = 0; N < 5e3; ++N)
    {
        std::cout << "\rmatrices multiplication... #" << N + 1 << std::flush;
        
        // generate random matrices
        int rowsA = rand() % 16 + 1;
        int colsArowsB = rand() % 16 + 1;
        int colsB = rand() % 16 + 1;

        std::vector<std::vector<int>> classicMatrixA = generateRandomMatrix<int>(rowsA, colsArowsB);
        CSRMatrixMock<int> sparseMatrixA = CSRMatrixMock<int>::fromVectors(classicMatrixA);

        std::vector<std::vector<int>> classicMatrixB = generateRandomMatrix<int>(colsArowsB, colsB);
        CSRMatrixMock<int> sparseMatrixB = CSRMatrixMock<int>::fromVectors(classicMatrixB);

        // calculate result manually
        std::vector<std::vector<int>> manualResult = multiplyMatrices<int>(classicMatrixA, classicMatrixB);

        // method
        assertEquals<SparseMatrix::CSRMatrix<int>, std::vector<std::vector<int>>>(
            sparseMatrixA.multiply(sparseMatrixB),
            manualResult,
            "Incorrect matrix-matrix multiplication result");

        // operator
        assertEquals<SparseMatrix::CSRMatrix<int>, std::vector<std::vector<int>>>(
            sparseMatrixA * sparseMatrixB,
            manualResult,
            "Incorrect matrix-matrix multiplication result (operator *)");
    }

    std::cout << " OK" << std::endl;
}