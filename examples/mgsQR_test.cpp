#include <iostream>
#include <vector>
#include <iomanip>

#include "mgsQR.hpp"

int main() {
    std::vector<std::vector<double>> a = {{1.0, 2.0, 3.0, 4.0}, {2.0, 4.0, 6.0, 5.0}, {3.0, 6.0, 9.0, 6.0}};
    // std::vector<double> rhs = {9.0, 15.0, 10.0, 12.0, 18.0};
    std::vector<double> x(3);

    // Container for the rectangular matrix R
    std::vector<std::vector<double>> r(a.size(), std::vector<double>(a.size(), 0.0));

    // Print the matrix A
    std::cout<< std::endl << "A (before factorization): " << std::endl;
    for (size_t i = 0; i < a[0].size(); i++) {
        for (size_t j = 0; j < a.size(); j++) {
            std::cout << std::setw(10) << a[j][i] << " ";
        }
        std::cout << std::endl;
    }

    // Print the matrix R
    std::cout<< std::endl << "R (before factorization): " << std::endl;
    for (size_t i = 0; i < r[0].size(); i++) {
        for (size_t j = 0; j < r.size(); j++) {
            std::cout << std::setw(10) << r[j][i] << " ";
        }
        std::cout << std::endl;
    }

    // perfrom mgs qr factorization
    gramSchmidt(a, r, 3, 3);

    // Solve the problem
    // householderQRSolve(a, rhs, x, 5, 3);

    // Print the matrix
    // Matrix A is converted to the Q matrix
    std::cout<< std::endl << "A (after factorization): " << std::endl;
    for (size_t i = 0; i < a[0].size(); i++) {
        for (size_t j = 0; j < a.size(); j++) {
            std::cout << std::setw(10) << a[j][i] << " ";
        }
        std::cout << std::endl;
    }

    // Print the matrix R
    std::cout<< std::endl << "R (after factorization): " << std::endl;
    for (size_t i = 0; i < r[0].size(); i++) {
        for (size_t j = 0; j < r.size(); j++) {
            std::cout << std::setw(10) << r[j][i] << " ";
        }
        std::cout << std::endl;
    }

    // Print the rhs
    /* std::cout << "rhs: " << std::endl;
    for (int i = 0; i < 3; i++) {
        std::cout << rhs[i] << "\n";
    }
    std::cout << std::endl; */

    // Print the solution
    /* std::cout << "x: " << std::endl;
    for (int i = 0; i < 3; i++) {
        std::cout << x[i] << "\n";
    }
    std::cout << std::endl; */
}

