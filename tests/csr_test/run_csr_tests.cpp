#include <iostream>

#include "cases/constructor.cpp"
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
        testConstructorFail1();
        testConstructorFail2();
        testConstructorFail3();
        testConstructorFail4();
        testConstructorFail5();
        testConstructorFail6();
        testConstructorFail7();
        testConstructorFail8();
        testConstructorFail9();
        testInternalStorage();
        testAdditionFail1();
        testAdditionFail2();
        testAdditionFail3();
        testAddition();
        // testMultiplicationFail1();
        // testMultiplicationFail2();
        testVectorMultiplication();
        // testMatrixMultiplication();
        testOutput();
        testSubtractionFail1();
        testSubtractionFail2();
        testSubtractionFail3();
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