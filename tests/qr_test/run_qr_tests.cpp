#include <iostream>

#include "cases/test_mgs.cpp"
#include "cases/test_householder.cpp"
#include "cases/test_eigen.cpp"

int main(int argc, char** argv)
{
    std::cout << "Running MGS QR tests..." << std::endl;
    testMgsQRSolve1();
    testMgsQRSolve2();
    testMgsQRSolve3();
    testMgsQRSolve4();
    testMgsQRSolve5();
    // testMgsQRSolve6();

    std::cout << "\nRunning Householder QR tests..." << std::endl;
    testHouseholderQRSolve1();
    testHouseholderQRSolve2();
    testHouseholderQRSolve3();
    testHouseholderQRSolve4();
    testHouseholderQRSolve5();
    testHouseholderQRSolve6();

    std::cout << "\nRunning Eigen QR tests..." << std::endl;
    testEigenQRSolve1();
    testEigenQRSolve2();
    testEigenQRSolve3();
    testEigenQRSolve4();
    testEigenQRSolve5();
    testEigenQRSolve6();
}