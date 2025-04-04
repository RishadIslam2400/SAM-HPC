#include "CSRMatrixMock.hpp"
#include "testlib.hpp"
#include "helpers.hpp"

void _subtractionFail1()
{
    SparseMatrix::CSRMatrix<int> a(3, 4), b(3, 5);
    a.subtract(b);
}

void testSubtractionFail1()
{
    std::cout << "subtract() fail #1..." << std::flush;
    assertException("InvalidDimensionsException", _subtractionFail1);
    std::cout << " OK" << std::endl;
}

void _subtractionFail2()
{
    SparseMatrix::CSRMatrix<int> a(3, 4), b(4, 4);
    a.subtract(b);
}

void testSubtractionFail2()
{
    std::cout << "subtract() fail #2..." << std::flush;
    assertException("InvalidDimensionsException", _subtractionFail2);
    std::cout << " OK" << std::endl;
}

void _subtractionFail3()
{
    SparseMatrix::CSRMatrix<int> a(3, 4), b(4, 5);
    a.subtract(b);
}

void testSubtractionFail3()
{
    std::cout << "subtract() fail #3..." << std::flush;
    assertException("InvalidDimensionsException", _subtractionFail3);
    std::cout << " OK" << std::endl;
}

void testSubtraction()
{
    for (int N = 0; N < 5e3; ++N)
    {
        std::cout << "\rmatrices subtraction...#" << N + 1 << std::flush;

        // generate random matrices
        int rows = rand() % 16 + 1;
        int cols = rand() % 16 + 1;

        std::vector<std::vector<int>> classicMatrixA = generateRandomMatrix<int>(rows, cols);
        CSRMatrixMock<int> sparseMatrixA = CSRMatrixMock<int>::fromVectors(classicMatrixA);

        std::vector<std::vector<int>> classicMatrixB = generateRandomMatrix<int>(rows, cols);
        CSRMatrixMock<int> sparseMatrixB = CSRMatrixMock<int>::fromVectors(classicMatrixB);

        // calculate results manually
        std::vector<std::vector<int>> manualResult = subtractMatrices<int>(classicMatrixA, classicMatrixB);

        // method
        assertEquals<SparseMatrix::CSRMatrix<int>, std::vector<std::vector<int>>>(
            sparseMatrixA.subtract(sparseMatrixB),
            manualResult,
            "incorrect matrices addition");

        // operator
        assertEquals<SparseMatrix::CSRMatrix<int>, std::vector<std::vector<int>>>(
            sparseMatrixA - sparseMatrixB,
            manualResult,
            "incorrect matrices addition (- operator)");
    }

    std::cout << " OK" << std::endl;
}