#include "CSRMatrixMock.hpp"

void testConstructorFail1()
{
    std::cout << "constructor fail #1..." << std::endl;
    SparseMatrix::CSRMatrix<int>(0);
    std::cout << "OK" << std::endl;
}