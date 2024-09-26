#include <iostream>
#include <vector>
#include <iomanip>

#include "householderQR.hpp"

int main() {
    std::vector<std::vector<double>> a = {{1.0, 3.0, 2.0, 1.0, 4.0}, {2.0, 2.0, 1.0, 3.0, 1.0}, {1.0, 1.0, 2.0, 2.0, 3.0}};
    std::vector<double> rhs = {9.0, 15.0, 10.0, 12.0, 18.0};
    std::vector<double> x(3);

    // Solve the problem
    householderQRSolve(a, rhs, x, 5, 3);

    // Print the matrix
    std::cout << "A: " << std::endl;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            std::cout << std::setw(10) << a[j][i] << " ";
        }
        std::cout << std::endl;
    }

    // Print the rhs
    std::cout << "rhs: " << std::endl;
    for (int i = 0; i < 3; i++) {
        std::cout << rhs[i] << "\n";
    }
    std::cout << std::endl;

    // Print the solution
    std::cout << "x: " << std::endl;
    for (int i = 0; i < 3; i++) {
        std::cout << x[i] << "\n";
    }
    std::cout << std::endl;
}

