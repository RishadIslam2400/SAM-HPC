#include "CSRMatrixMock.hpp"
#include "testlib.hpp"
#include "helpers.hpp"

void testSubtraction()
{
    for (int N = 0; N < 5e3; ++N)
    {
        std::cout << "\rmatrices subtraction...#" << N + 1 << std::flush;

        // generate random matrices
        size_t rows = rand() % 16 + 1;
        size_t cols = rand() % 16 + 1;

        std::vector<std::vector<int>> classicMatrixA = generateRandomMatrix<int>(rows, cols);
        CSRMatrixMock<int> sparseMatrixA(classicMatrixA);

        std::vector<std::vector<int>> classicMatrixB = generateRandomMatrix<int>(rows, cols);
        CSRMatrixMock<int> sparseMatrixB(classicMatrixB);

        // calculate results manually
        std::vector<std::vector<int>> manualResult = subtractMatrices<int>(classicMatrixA, classicMatrixB);

        // method
        assertEquals<CSRMatrix<int>, std::vector<std::vector<int>>>(
            *sparseMatrixA.subtract(sparseMatrixB),
            manualResult,
            "incorrect matrices subtraction");

        // operator
        assertEquals<CSRMatrix<int>, std::vector<std::vector<int>>>(
            *(sparseMatrixA - sparseMatrixB),
            manualResult,
            "incorrect matrices subtraction (+ operator)");
    }

    std::cout << " OK" << std::endl;
}