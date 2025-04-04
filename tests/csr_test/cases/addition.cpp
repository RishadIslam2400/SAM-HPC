#include "CSRMatrixMock.hpp"
#include "testlib.hpp"
#include "helpers.hpp"

void _additionFail1()
{
    SparseMatrix::CSRMatrix<int> a(3, 4), b(3, 5);
    a.add(b);
}

void testAdditionFail1()
{
    std::cout << "add() fail #1..." << std::flush;
    assertException("InvalidDimensionsException", _additionFail1);
    std::cout << " OK" << std::endl;
}

void _additionFail2()
{
    SparseMatrix::CSRMatrix<int> a(3, 4), b(4, 4);
    a.add(b);
}

void testAdditionFail2()
{
    std::cout << "add() fail #2..." << std::flush;
    assertException("InvalidDimensionsException", _additionFail2);
    std::cout << " OK" << std::endl;
}

void _additionFail3()
{
    SparseMatrix::CSRMatrix<int> a(3, 4), b(4, 5);
    a.add(b);
}

void testAdditionFail3()
{
    std::cout << "add() fail #3..." << std::flush;
    assertException("InvalidDimensionsException", _additionFail3);
    std::cout << " OK" << std::endl;
}

void testAddition()
{
    for (int N = 0; N < 5e3; ++N)
    {
        std::cout << "\rmatrices addition...#" << N + 1 << std::flush;

        // generate random matrices
        int rows = rand() % 16 + 1;
        int cols = rand() % 16 + 1;

        std::vector<std::vector<int>> classicMatrixA = generateRandomMatrix<int>(rows, cols);
        CSRMatrixMock<int> sparseMatrixA = CSRMatrixMock<int>::fromVectors(classicMatrixA);

        std::vector<std::vector<int>> classicMatrixB = generateRandomMatrix<int>(rows, cols);
        CSRMatrixMock<int> sparseMatrixB = CSRMatrixMock<int>::fromVectors(classicMatrixB);

        // calculate results manually
        std::vector<std::vector<int>> manualResult = addMatrices<int>(classicMatrixA, classicMatrixB);

        // method
        assertEquals<SparseMatrix::CSRMatrix<int>, std::vector<std::vector<int>>>(
            sparseMatrixA.add(sparseMatrixB),
            manualResult,
            "incorrect matrices addition");

        // operator
        assertEquals<SparseMatrix::CSRMatrix<int>, std::vector<std::vector<int>>>(
            sparseMatrixA + sparseMatrixB,
            manualResult,
            "incorrect matrices addition (+ operator)");
    }

    std::cout << " OK" << std::endl;
}