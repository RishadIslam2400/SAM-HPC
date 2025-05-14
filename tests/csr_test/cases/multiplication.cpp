#include "CSRMatrixMock.hpp"
#include "testlib.hpp"

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