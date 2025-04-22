#include <iostream>

#include "cases/test_mgs.cpp"

int main(int argc, char** argv)
{
    std::cout << "Running QR tests..." << std::endl;
    testMgsQRSolve1();
    testMgsQRSolve2();
    testMgsQRSolve3();
    testMgsQRSolve4();
    testMgsQRSolve5();
    testMgsQRSolve6();
}