#include <iostream>

#include "cases/addition.cpp"
#include "cases/csr_format.cpp"
#include "cases/multiplication.cpp"
#include "cases/output.cpp"
#include "cases/subtraction.cpp"
#include "cases/linear_algebra.cpp"

int main(int argc, char** argv)
{
    std::cout << "Running CSR tests..." << std::endl;
    srand(time(NULL));

    try
    {
        testInternalStorage();
        testAddition();
        testVectorMultiplication();
        testMatrixMultiplication();
        testOutput();
        testSubtraction();
        testSortRows();
        testTranspose();
        testDiagonal();
    }
    catch(const FailureException& e)
    {
        std::cout << "- Fail: '" << e.what() << "'" << std::endl;
    }

    return 0;
}