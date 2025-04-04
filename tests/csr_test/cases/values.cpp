#include "CSRMatrixMock.hpp"
#include "testlib.hpp"

void _getFail()
{
    SparseMatrix::CSRMatrix<int> m(1, 1);
    m.get(1, 0);
}

void testGetFail()
{
    std::cout << "get() fail..." << std::flush;
    assertException("InvalidCoordinatesException", _getFail);
    std::cout << "OK" << std::endl;
}

void _setFail()
{
    SparseMatrix::CSRMatrix<int> m(3, 4);
    m.set(-1, 3, 0);
}

void testSetFail()
{
    std::cout << "set() fail..." << std::flush;
    assertException("InvalidCoordinatesException", _setFail);
    std::cout << "OK" << std::endl;
}

void testGettersAndSettters()
{
    std::cout << "getters/settters..." << std::flush;

    SparseMatrix::CSRMatrix<int> m(3);
    for (int i = 0; i < 3; ++i) 
    {
        for (int j = 0; j < 3; ++j) 
        {
            assertEquals<int>(0, m.get(i, j));
        }
    }

    m.set(-4, 0, 2);
    assertEquals<int>(-4, m.get(0, 2));

    std::cout << "OK" << std::endl;
}