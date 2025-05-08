#include "CSRMatrixMock.hpp"
#include "testlib.hpp"
#include "helpers.hpp"

void _subtractionFail1()
{
    // Create two matrices with different dimensions
    /*
            [0 2 0 1]                                    [6 7 0]
        A = [1 3 0 0]                                B = [0 1 0]
            [0 7 4 0]                                    [0 0 2]
        row pointers:   [0, 2, 4, 6]                     [7 3 2]
        column indices: [1, 3, 0, 1, 1, 2]    row pointers:   [0, 2, 3, 4, 8]
        values:         [2, 1, 1, 3, 7, 4]    column indices: [0, 1, 1, 2, 0, 1, 2]
                                              values:         [6, 7, 1, 2, 7, 3, 2]
    */
    std::vector<size_t> rowPointersA = {0, 2, 4, 6};
    std::vector<size_t> colIndicesA = {1, 3, 0, 1, 1, 2};
    std::vector<int> valsA = {2, 1, 1, 3, 7, 4};
    CSRMatrix<int> A(3, 4, valsA, rowPointersA, colIndicesA);

    std::vector<size_t> rowPointersB = {0, 2, 3, 4, 8};
    std::vector<size_t> colIndicesB = {0, 1, 1, 2, 0, 1, 2};
    std::vector<int> valsB = {6, 7, 1, 2, 7, 3, 2};
    CSRMatrix<int> B(4, 3, valsB, rowPointersB, colIndicesB);

    A.subtract(B);
}

void testSubtractionFail1()
{
    std::cout << "add() fail #1..." << std::flush;
    assertException("InvalidDimensionsException", _subtractionFail1);
    std::cout << " OK" << std::endl;
}

void _subtractionFail2()
{
    // Create two matrices with different dimensions
    /*
            [0 2 0 1]                                    [6 7 0 7]
        A = [1 3 0 0]                                B = [0 1 0 0]
            [0 7 4 0]                                    [0 0 2 0]
        row pointers:   [0, 2, 4, 6]                     [7 3 2 2]
        column indices: [1, 3, 0, 1, 1, 2]    row pointers:   [0, 3, 4, 5, 9]
        values:         [2, 1, 1, 3, 7, 4]    column indices: [0, 1, 2, 1, 2, 0, 1, 2, 3]
                                              values:         [6, 7, 7, 1, 2, 7, 3, 2, 2]
    */
    std::vector<size_t> rowPointersA = {0, 2, 4, 6};
    std::vector<size_t> colIndicesA = {1, 3, 0, 1, 1, 2};
    std::vector<int> valsA = {2, 1, 1, 3, 7, 4};
    CSRMatrix<int> A(3, 4, valsA, rowPointersA, colIndicesA);

    std::vector<size_t> rowPointersB = {0, 3, 4, 5, 9};
    std::vector<size_t> colIndicesB = {0, 1, 2, 1, 2, 0, 1, 2, 3};
    std::vector<int> valsB = {6, 7, 7, 1, 2, 7, 3, 2, 2};
    CSRMatrix<int> B(4, 4, valsB, rowPointersB, colIndicesB);

    A.subtract(B);
}

void testSubtractionFail2()
{
    std::cout << "add() fail #2..." << std::flush;
    assertException("InvalidDimensionsException", _subtractionFail2);
    std::cout << " OK" << std::endl;
}

void _subtractionFail3()
{
    // Create two matrices with different dimensions
    /*
            [0 2 0 1]                                    [6 7 0 7 0]
        A = [1 3 0 0]                                B = [0 1 0 0 0]
            [0 7 4 0]                                    [0 0 2 0 0]
        row pointers:   [0, 2, 4, 6]                     [7 3 2 2 0]
        column indices: [1, 3, 0, 1, 1, 2]    row pointers:   [0, 3, 4, 5, 9]
        values:         [2, 1, 1, 3, 7, 4]    column indices: [0, 1, 2, 1, 2, 0, 1, 2, 3]
                                              values:         [6, 7, 7, 1, 2, 7, 3, 2, 2]
    */
    std::vector<size_t> rowPointersA = {0, 2, 4, 6};
    std::vector<size_t> colIndicesA = {1, 3, 0, 1, 1, 2};
    std::vector<int> valsA = {2, 1, 1, 3, 7, 4};
    CSRMatrix<int> A(3, 4, valsA, rowPointersA, colIndicesA);

    std::vector<size_t> rowPointersB = {0, 3, 4, 5, 9};
    std::vector<size_t> colIndicesB = {0, 1, 2, 1, 2, 0, 1, 2, 3};
    std::vector<int> valsB = {6, 7, 7, 1, 2, 7, 3, 2, 2};
    CSRMatrix<int> B(4, 5, valsB, rowPointersB, colIndicesB);

    A.subtract(B);
}

void testSubtractionFail3()
{
    std::cout << "add() fail #3..." << std::flush;
    assertException("InvalidDimensionsException", _subtractionFail3);
    std::cout << " OK" << std::endl;
}

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