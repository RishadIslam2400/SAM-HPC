#include <iostream>
#include <chrono>

#include "cases/constructor.cpp"

int main(int argc, char** argv)
{
    std::cout << "Running CSR tests..." << std::endl;
    srand(time(NULL));

    testConstructorFail1();
    return 0;
}