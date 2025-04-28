#include "householderQR.hpp"
#include "mgsQR.hpp"
#include "eigenQRSolve.hpp"

#include <iostream>
#include <algorithm>
#include <random>
#include <chrono>

std::vector<std::vector<double>> generateRandomMatrix(const size_t rows, const size_t cols, const double min = -100.0, const double max = 100.0)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(min, max);
    std::vector<std::vector<double>> matrix(cols, std::vector<double>(rows));

    for (size_t i = 0; i < cols; ++i)
    {
        for (size_t j = 0; j < rows; ++j)
        {
            matrix[i][j] = dis(gen);
        }
    }

    return matrix;
}

std::vector<double> generateRandomVector(const size_t size, const double min = -100.0, const double max = 100.0)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(min, max);
    std::vector<double> vector(size);

    for (size_t i = 0; i < size; ++i)
    {
        vector[i] = dis(gen);
    }

    return vector;
}

int main(int argc, char **argv)
{
    srand(time(0));
    const size_t maxIter = 50;

    std::vector<std::pair<size_t, size_t>> tallSkinnyMatrixDimensions = {
        {100, 25},
        {200, 50},
        {375, 53},
        {500, 100},
        {1000, 300}
    };

    std::vector<std::pair<size_t, size_t>> squareMatrixDimensions = {
        {100, 100},
        {200, 200},
        {375, 375},
        {500, 500},
        {1000, 1000}
    };

    std::vector<std::pair<size_t, size_t>> wideMatrixDimensions = {
        {200, 150},
        {375, 200},
        {500, 300},
        {1000, 700},
    };

    std::vector<std::pair<size_t, size_t>> smallMatrixDimensions = {
        {20, 5},
        {50, 10},
        {40, 25},
        {80, 20},
        {60, 30}
    };

    std::cout << "Testing tall skinny matrices:" << std::endl;
    for (size_t i = 0; i < tallSkinnyMatrixDimensions.size(); ++i)
    {
        std::vector<std::vector<double>> testMatrix = generateRandomMatrix(tallSkinnyMatrixDimensions[i].first, tallSkinnyMatrixDimensions[i].second);
        std::vector<double> testRHS = generateRandomVector(tallSkinnyMatrixDimensions[i].first);
        std::vector<double> x(tallSkinnyMatrixDimensions[i].second);
        std::cout << "Matrix size: " << tallSkinnyMatrixDimensions[i].first << " x " << tallSkinnyMatrixDimensions[i].second << std::endl;


        std::chrono::high_resolution_clock::time_point mgsStart = std::chrono::high_resolution_clock::now();
        for (size_t iter = 0; iter < maxIter; ++iter)
        {
            mgsQRSolve(testMatrix, testRHS, x, tallSkinnyMatrixDimensions[i].first, tallSkinnyMatrixDimensions[i].second);
        }
        std::chrono::high_resolution_clock::time_point mgsEnd = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> mgsElapsed = mgsEnd - mgsStart;
        std::cout << "MGS QR Solve: " << mgsElapsed.count() / maxIter << "s" << std::endl;

        std::chrono::high_resolution_clock::time_point householderStart = std::chrono::high_resolution_clock::now();
        for (size_t iter = 0; iter < maxIter; ++iter)
        {
            householderQRSolve(testMatrix, testRHS, x, tallSkinnyMatrixDimensions[i].first, tallSkinnyMatrixDimensions[i].second);
        }
        std::chrono::high_resolution_clock::time_point householderEnd = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> householderElapsed = householderEnd - householderStart;
        std::cout << "Householder QR Solve: " << householderElapsed.count() / maxIter << "s" << std::endl;

        std::chrono::high_resolution_clock::time_point eigenStart = std::chrono::high_resolution_clock::now();
        for (size_t iter = 0; iter < maxIter; ++iter)
        {
            eigenQRSolve(testMatrix, testRHS, x, tallSkinnyMatrixDimensions[i].first, tallSkinnyMatrixDimensions[i].second);
        }
        std::chrono::high_resolution_clock::time_point eigenEnd = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> eigenElapsed = eigenEnd - eigenStart;
        std::cout << "Eigen QR Solve: " << eigenElapsed.count() / maxIter << "s" << std::endl;
    }

    std::cout << "Testing square matrices:" << std::endl;
    for (size_t i = 0; i < squareMatrixDimensions.size(); ++i)
    {
        std::vector<std::vector<double>> testMatrix = generateRandomMatrix(squareMatrixDimensions[i].first, squareMatrixDimensions[i].second);
        std::vector<double> testRHS = generateRandomVector(squareMatrixDimensions[i].first);
        std::vector<double> x(squareMatrixDimensions[i].second);
        std::cout << "Matrix size: " << squareMatrixDimensions[i].first << " x " << squareMatrixDimensions[i].second << std::endl;

        std::chrono::high_resolution_clock::time_point mgsStart = std::chrono::high_resolution_clock::now();
        for (size_t iter = 0; iter < maxIter; ++iter)
        {
            mgsQRSolve(testMatrix, testRHS, x, squareMatrixDimensions[i].first, squareMatrixDimensions[i].second);
        }
        std::chrono::high_resolution_clock::time_point mgsEnd = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> mgsElapsed = mgsEnd - mgsStart;
        std::cout << "MGS QR Solve: " << mgsElapsed.count() / maxIter << "s" << std::endl;

        std::chrono::high_resolution_clock::time_point householderStart = std::chrono::high_resolution_clock::now();
        for (size_t iter = 0; iter < maxIter; ++iter)
        {
            householderQRSolve(testMatrix, testRHS, x, squareMatrixDimensions[i].first, squareMatrixDimensions[i].second);
        }
        std::chrono::high_resolution_clock::time_point householderEnd = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> householderElapsed = householderEnd - householderStart;
        std::cout << "Householder QR Solve: " << householderElapsed.count() / maxIter << "s" << std::endl;

        std::chrono::high_resolution_clock::time_point eigenStart = std::chrono::high_resolution_clock::now();
        for (size_t iter = 0; iter < maxIter; ++iter)
        {
            eigenQRSolve(testMatrix, testRHS, x, squareMatrixDimensions[i].first, squareMatrixDimensions[i].second);
        }
        std::chrono::high_resolution_clock::time_point eigenEnd = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> eigenElapsed = eigenEnd - eigenStart;
        std::cout << "Eigen QR Solve: " << eigenElapsed.count() / maxIter << "s" << std::endl;
    }

    std::cout << "Testing wide matrices:" << std::endl;
    for (size_t i = 0; i < wideMatrixDimensions.size(); ++i)
    {
        std::vector<std::vector<double>> testMatrix = generateRandomMatrix(wideMatrixDimensions[i].first, wideMatrixDimensions[i].second);
        std::vector<double> testRHS = generateRandomVector(wideMatrixDimensions[i].first);
        std::vector<double> x(wideMatrixDimensions[i].second);
        std::cout << "Matrix size: " << wideMatrixDimensions[i].first << " x " << wideMatrixDimensions[i].second << std::endl;

        std::chrono::high_resolution_clock::time_point mgsStart = std::chrono::high_resolution_clock::now();
        for (size_t iter = 0; iter < maxIter; ++iter)
        {
            mgsQRSolve(testMatrix, testRHS, x, wideMatrixDimensions[i].first, wideMatrixDimensions[i].second);
        }
        std::chrono::high_resolution_clock::time_point mgsEnd = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> mgsElapsed = mgsEnd - mgsStart;
        std::cout << "MGS QR Solve: " << mgsElapsed.count() / maxIter << "s" << std::endl;

        std::chrono::high_resolution_clock::time_point householderStart = std::chrono::high_resolution_clock::now();
        for (size_t iter = 0; iter < maxIter; ++iter)
        {
            householderQRSolve(testMatrix, testRHS, x, wideMatrixDimensions[i].first, wideMatrixDimensions[i].second);
        }
        std::chrono::high_resolution_clock::time_point householderEnd = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> householderElapsed = householderEnd - householderStart;
        std::cout << "Householder QR Solve: " << householderElapsed.count() / maxIter << "s" << std::endl;

        std::chrono::high_resolution_clock::time_point eigenStart = std::chrono::high_resolution_clock::now();
        for (size_t iter = 0; iter < maxIter; ++iter)
        {
            eigenQRSolve(testMatrix, testRHS, x, wideMatrixDimensions[i].first, wideMatrixDimensions[i].second);
        }
        std::chrono::high_resolution_clock::time_point eigenEnd = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> eigenElapsed = eigenEnd - eigenStart;
        std::cout << "Eigen QR Solve: " << eigenElapsed.count() / maxIter << "s" << std::endl;
    }

    std::cout << "Testing small matrices:" << std::endl;
    for (size_t i = 0; i < smallMatrixDimensions.size(); ++i)
    {
        std::vector<std::vector<double>> testMatrix = generateRandomMatrix(smallMatrixDimensions[i].first, smallMatrixDimensions[i].second);
        std::vector<double> testRHS = generateRandomVector(smallMatrixDimensions[i].first);
        std::vector<double> x(smallMatrixDimensions[i].second);
        std::cout << "Matrix size: " << smallMatrixDimensions[i].first << " x " << smallMatrixDimensions[i].second << std::endl;

        std::chrono::high_resolution_clock::time_point mgsStart = std::chrono::high_resolution_clock::now();
        for (size_t iter = 0; iter < maxIter; ++iter)
        {
            mgsQRSolve(testMatrix, testRHS, x, smallMatrixDimensions[i].first, smallMatrixDimensions[i].second);
        }
        std::chrono::high_resolution_clock::time_point mgsEnd = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> mgsElapsed = mgsEnd - mgsStart;
        std::cout << "MGS QR Solve: " << mgsElapsed.count() / maxIter << "s" << std::endl;

        std::chrono::high_resolution_clock::time_point householderStart = std::chrono::high_resolution_clock::now();
        for (size_t iter = 0; iter < maxIter; ++iter)
        {
            householderQRSolve(testMatrix, testRHS, x, smallMatrixDimensions[i].first, smallMatrixDimensions[i].second);
        }
        std::chrono::high_resolution_clock::time_point householderEnd = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> householderElapsed = householderEnd - householderStart;
        std::cout << "Householder QR Solve: " << householderElapsed.count() / maxIter << "s" << std::endl;

        std::chrono::high_resolution_clock::time_point eigenStart = std::chrono::high_resolution_clock::now();
        for (size_t iter = 0; iter < maxIter; ++iter)
        {
            eigenQRSolve(testMatrix, testRHS, x, smallMatrixDimensions[i].first, smallMatrixDimensions[i].second);
        }
        std::chrono::high_resolution_clock::time_point eigenEnd = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> eigenElapsed = eigenEnd - eigenStart;
        std::cout << "Eigen QR Solve: " << eigenElapsed.count() / maxIter << "s" << std::endl;
    }
}