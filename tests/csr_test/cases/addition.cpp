#include "CSRMatrixMock.hpp"
#include "testlib.hpp"
#include "helpers.hpp"

void testAddition()
{
    for (int N = 0; N < 5e3; ++N)
    {
        std::cout << "\rmatrices addition...#" << N + 1 << std::flush;

        // generate random matrices
        size_t rows = rand() % 16 + 1;
        size_t cols = rand() % 16 + 1;

        std::vector<std::vector<int>> classicMatrixA = generateRandomMatrix<int>(rows, cols);
        CSRMatrixMock<int> sparseMatrixA(classicMatrixA);

        std::vector<std::vector<int>> classicMatrixB = generateRandomMatrix<int>(rows, cols);
        CSRMatrixMock<int> sparseMatrixB(classicMatrixB);

        // calculate results manually
        std::vector<std::vector<int>> manualResult = addMatrices<int>(classicMatrixA, classicMatrixB);

        // method
        assertEquals<CSRMatrix<int>, std::vector<std::vector<int>>>(
            *(sparseMatrixA.add(sparseMatrixB)),
            manualResult,
            "incorrect matrices addition");

        // operator
        assertEquals<CSRMatrix<int>, std::vector<std::vector<int>>>(
            *(sparseMatrixA + sparseMatrixB),
            manualResult,
            "incorrect matrices addition (+ operator)");
    }

    std::cout << " OK" << std::endl;
}