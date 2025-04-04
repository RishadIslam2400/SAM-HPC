#include "CSRMatrixMock.hpp"
#include "testlib.hpp"

void _constructorFail1()
{
    SparseMatrix::CSRMatrix<int>(0);
}

void testConstructorFail1()
{
    std::cout << "constructor fail #1..." << std::flush;
    assertException("InvalidDimensionsException", _constructorFail1);
    std::cout << "OK" << std::endl;
}

void _constructorFail2()
{
    SparseMatrix::CSRMatrix<int>(0, 1);
}

void testConstructorFail2()
{
    std::cout << "constructor fail #2..." << std::flush;
    assertException("InvalidDimensionsException", _constructorFail2);
    std::cout << " OK" << std::endl;
}

void _constructorFail3()
{
    SparseMatrix::CSRMatrix<int>(1, 0);
}

void testConstructorFail3()
{
    std::cout << "constructor fail #3..." << std::flush;
    assertException("InvalidDimensionsException", _constructorFail3);
    std::cout << " OK" << std::endl;
}

void _constructorFail4()
{
    SparseMatrix::CSRMatrix<int>(0, 0);
}

void testConstructorFail4()
{
    std::cout << "constructor fail #4..." << std::flush;
    assertException("InvalidDimensionsException", _constructorFail4);
    std::cout << " OK" << std::endl;
}