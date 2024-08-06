#ifndef HELPER_FUNCTIONS_H
#define HELPER_FUNCTIONS_H

#include <cusolverDn.h>
#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <iostream>
#include <iomanip>

// Utility function to check CUDA errors
inline void checkCuda(cudaError_t result, const char* msg) {
    if (result != cudaSuccess) {
        std::cerr << "CUDA Error: " << msg << " - " << cudaGetErrorString(result) << std::endl;
        exit(EXIT_FAILURE);
    }
}

// Utility function to check cuSOLVER errors
inline void checkCusolver(cusolverStatus_t result, const char* msg) {
    if (result != CUSOLVER_STATUS_SUCCESS) {
        std::cerr << "cuSOLVER Error: " << msg << std::endl;
        exit(EXIT_FAILURE);
    }
}

// Utility function to check cuBLAS errors
inline void checkCublas(cublasStatus_t result, const char* msg) {
    if (result != CUBLAS_STATUS_SUCCESS) {
        std::cerr << "cuBLAS Error: " << msg << std::endl;
        exit(EXIT_FAILURE);
    }
}

// Function to print a matrix
inline void printMatrix(const char* name, double* matrix, int rows, int cols) {
    std::cout << name << " Matrix:" << std::endl;
    for (int j = 0; j < rows; ++j) {
        for (int i = 0; i < cols; ++i) {
            double value = matrix[j * cols + i];
            std::cout << std::setw(5) << value << " ";
        }
        std::cout << std::endl;
    }
}

// Function to print a vector
inline void printVector(const char* name, double* vector, int size) {
    std::cout << name << " Vector:" << std::endl;
    for (int i = 0; i < size; ++i) {
        std::cout << std::setw(5) << vector[i] << " ";
    }
    std::cout << std::endl;
}

// Function to print the contents of a device array
inline void printDeviceArray(const char* name, double* d_array, int size) {
    double* h_array = (double*)malloc(size * sizeof(double));
    checkCuda(cudaMemcpy(h_array, d_array, size * sizeof(double), cudaMemcpyDeviceToHost), "Failed to copy device array to host");
    std::cout << name << ":" << std::endl;
    for (int i = 0; i < size; ++i) {
        std::cout << h_array[i] << " ";
    }
    std::cout << std::endl;
    free(h_array);
}

#endif // HELPER_FUNCTIONS_H
