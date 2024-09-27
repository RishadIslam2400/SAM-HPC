#include <iostream>
#include <vector>
#include <iomanip>

#include "mgsQR.hpp"

int main() {
    std::vector<std::vector<double>> a = {{1.0, 1.0, 1.0, 1.0}, {1.0, 1.00001, 1.0, 1.0}, {1.0, 1.0, 1.00001, 1.0}, {1.0, 1.0, 1.0, 1.00001}};
    std::vector<double> rhs = {4.0, 4.00001, 4.00001, 4.00001};
    std::vector<double> x(a.size(), 0.0);

    // Print the matrix A
    std::cout<< std::endl << "A (before factorization): " << std::endl;
    for (size_t i = 0; i < a[0].size(); i++) {
        for (size_t j = 0; j < a.size(); j++) {
            std::cout << std::setw(10) << a[j][i] << " ";
        }
        std::cout << std::endl;
    }

    // perfrom mgs qr solve
    mgsQRSolve(a, rhs, x, a[0].size(), a.size());

    // Print the matrix
    // Matrix A is converted to the Q matrix
    std::cout<< std::endl << "A (after factorization): " << std::endl;
    for (size_t i = 0; i < a[0].size(); i++) {
        for (size_t j = 0; j < a.size(); j++) {
            std::cout << std::setw(10) << a[j][i] << " ";
        }
        std::cout << std::endl;
    }

    // Print the rhs
    std::cout << "rhs: " << std::endl;
    for (size_t i = 0; i < a[0].size(); i++) {
        std::cout << rhs[i] << "\n";
    }
    std::cout << std::endl;

    // Print the solution
    std::cout << "x: " << std::endl;
    for (size_t i = 0; i < a.size(); i++) {
        std::cout << x[i] << "\n";
    }
    std::cout << std::endl;
}

