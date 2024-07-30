#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cuda_runtime.h>
#include <stdio.h>
#include <thrust/device_vector.h>
#include <thrust/extrema.h>

__global__ void getSkAndRkMaxSize(int *S_row_indices, int *S_col_ptr, int *A_col_ptr, int *Sk_max_sizes, int *Rk_max_sizes, int n);
__global__ void extractSubMatrix(double *A_values, int *A_row_indices, int *A_col_ptr, 
                                 int *S_row_indices, int *S_col_ptr, 
                                 int *Sk_max_sizes, int *Rk_actual_sizes, 
                                 int *rk_idx, int *submatrix, int max_Sk_size, int max_Rk_size, int n);



const int N = 100; // Matrix size (example with N=100 for demonstration)



void readCSCFromFile(const std::string& filename,
                     std::vector<double>& values,
                     std::vector<int>& rowIndices,
                     std::vector<int>& colPointers) {
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Error opening file for reading" << std::endl;
        return;
    }

    int numRows, numCols, numNonZeros;
    infile >> numRows >> numCols >> numNonZeros;

    values.resize(numNonZeros);
    rowIndices.resize(numNonZeros);
    colPointers.resize(numCols + 1);

    for (int i = 0; i < numNonZeros; ++i) {
        infile >> values[i];
    }
    for (int i = 0; i < numNonZeros; ++i) {
        infile >> rowIndices[i];
    }
    for (int i = 0; i < numCols + 1; ++i) {
        infile >> colPointers[i];
    }

    infile.close();
}

void printCSCArrays(const std::vector<double>& values,
                    const std::vector<int>& rowIndices,
                    const std::vector<int>& colPointers) {
    std::cout << "Values Array: ";
    for (const auto& val : values) {
        std::cout << val << " ";
    }
    std::cout << std::endl;

    std::cout << "Row Indices Array: ";
    for (const auto& row : rowIndices) {
        std::cout << row << " ";
    }
    std::cout << std::endl;

    std::cout << "Column Pointers Array: ";
    for (const auto& col : colPointers) {
        std::cout << col << " ";
    }
    std::cout << std::endl;
}



int main() {
    std::vector<double> values_A, values_S;
    std::vector<int> rowIndices_A, rowIndices_S;
    std::vector<int> colPointers_A, colPointers_S;
    
    // Generate filename based on the matrix size
    std::string filename = "matrix_" + std::to_string(N) + ".csc";

    // Read the matrix from a file for both A and S
    readCSCFromFile(filename, values_A, rowIndices_A, colPointers_A);
    readCSCFromFile(filename, values_S, rowIndices_S, colPointers_S);

    // Allocate device memory
    int *d_rowIndices_A, *d_colPointers_A;
    int *d_rowIndices_S, *d_colPointers_S;
    int *d_Sk_max_sizes, *d_Rk_max_sizes;
    int *d_Rk_actual_sizes;

    cudaMalloc(&d_rowIndices_A, rowIndices_A.size() * sizeof(int));
    cudaMalloc(&d_colPointers_A, colPointers_A.size() * sizeof(int));
    cudaMalloc(&d_rowIndices_S, rowIndices_S.size() * sizeof(int));
    cudaMalloc(&d_colPointers_S, colPointers_S.size() * sizeof(int));
    cudaMalloc(&d_Sk_max_sizes, N * sizeof(int));
    cudaMalloc(&d_Rk_max_sizes, N * sizeof(int));
    cudaMalloc(&d_Rk_actual_sizes, N * sizeof(int));

    cudaMemcpy(d_rowIndices_A, rowIndices_A.data(), rowIndices_A.size() * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_colPointers_A, colPointers_A.data(), colPointers_A.size() * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_rowIndices_S, rowIndices_S.data(), rowIndices_S.size() * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_colPointers_S, colPointers_S.data(), colPointers_S.size() * sizeof(int), cudaMemcpyHostToDevice);


    // Launch the getSkAndRkMaxSize kernel
    int blockSize = 256;
    int numBlocks = (N + blockSize - 1) / blockSize;
    getSkAndRkMaxSize<<<numBlocks, blockSize>>>(d_rowIndices_S, d_colPointers_S, d_colPointers_A, d_Sk_max_sizes, d_Rk_max_sizes, N);

    // Use Thrust to find the maximum values in d_Sk_max_sizes and d_Rk_max_sizes
    thrust::device_ptr<int> d_Sk_max_sizes_ptr(d_Sk_max_sizes);
    thrust::device_ptr<int> d_Rk_max_sizes_ptr(d_Rk_max_sizes);

    int Sk_max_size = *(thrust::max_element(d_Sk_max_sizes_ptr, d_Sk_max_sizes_ptr + N));
    int Rk_max_size = *(thrust::max_element(d_Rk_max_sizes_ptr, d_Rk_max_sizes_ptr + N));

    // Allocate memory for rk_idx and submatrix
    int *d_rk_idx, *d_submatrix;
    cudaMalloc(&d_rk_idx, N * Sk_max_size * sizeof(int));
    cudaMalloc(&d_submatrix, N * Rk_max_size * Sk_max_size * sizeof(int));

    // Launch the extractSubMatrix kernel
    extractSubMatrix<<<numBlocks, blockSize>>>(values_A.data(), d_rowIndices_A, d_colPointers_A,
                                               d_rowIndices_S, d_colPointers_S,
                                               d_Sk_max_sizes, d_Rk_actual_sizes,
                                               d_rk_idx, d_submatrix, Sk_max_size, Rk_max_size, N);

    // Copy rk_idx, d_Rk_actual_sizes, and submatrix back to host
    std::vector<int> h_rk_idx(N * Sk_max_size);
    std::vector<int> h_Rk_actual_sizes(N);
    std::vector<int> h_submatrix(N * Rk_max_size * Sk_max_size);

    cudaMemcpy(h_rk_idx.data(), d_rk_idx, N * Sk_max_size * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_Rk_actual_sizes.data(), d_Rk_actual_sizes, N * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_submatrix.data(), d_submatrix, N * Rk_max_size * Sk_max_size * sizeof(int), cudaMemcpyDeviceToHost);


    // Print the rk_idx results
    for (int i = 0; i < N; ++i) {
        std::cout << "rk_idx for column " << i << ": ";
        for (int j = 0; j < Rk_max_size; ++j) {
            if (h_rk_idx[i * Rk_max_size + j] == -1) break; // End of valid entries
            std::cout << h_rk_idx[i * Rk_max_size + j] << " ";
        }
        std::cout << std::endl;
    }
    
    // Print the Rk_actual_sizes results
    // std::cout << "Rk_actual_sizes:" << std::endl;
    // for (int i = 0; i < N; ++i) {
    //     std::cout << "Column " << i << ": " << h_Rk_actual_sizes[i] << std::endl;
    // }

    // Print the submatrix results
    std::cout << "Submatrix:" << std::endl;
    for (int i = 0; i < N; ++i) {
        std::cout << "Submatrix for column " << i << ": " << std::endl;
        for (int j = 0; j < Rk_max_size; ++j) {
            for (int k = 0; k < Sk_max_size; ++k) {
                int value = h_submatrix[i * Rk_max_size * Sk_max_size + j * Sk_max_size + k];
                if (value == -1) break; // End of valid entries
                std::cout << value << " ";
            }
            std::cout << std::endl;
        }
    }

    // Free device memory
    cudaFree(d_rowIndices_A);
    cudaFree(d_colPointers_A);
    cudaFree(d_rowIndices_S);
    cudaFree(d_colPointers_S);
    cudaFree(d_Sk_max_sizes);
    cudaFree(d_Rk_max_sizes);
    cudaFree(d_Rk_actual_sizes);
    cudaFree(d_rk_idx);
    cudaFree(d_submatrix);

    return 0;
}

__global__ void getSkAndRkMaxSize(int *S_row_indices, int *S_col_ptr, int *A_col_ptr, int *Sk_max_sizes, int *Rk_max_sizes, int n) {
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    if (col < n) {
        int Sk_max_size = 0;
        int Rk_max_size = 0;
        int start_S = S_col_ptr[col];
        int end_S = S_col_ptr[col + 1];
        Sk_max_size = end_S - start_S;  // Number of non-zero entries in column `col` of S

        for (int i = start_S; i < end_S; ++i) {
            int row = S_row_indices[i];
            int start_A = A_col_ptr[row];
            int end_A = A_col_ptr[row + 1];
            Rk_max_size += (end_A - start_A);
        }
        
        Sk_max_sizes[col] = Sk_max_size;
        Rk_max_sizes[col] = Rk_max_size;
    }
}

__global__ void extractSubMatrix(double *A_values, int *A_row_indices, int *A_col_ptr, 
                                 int *S_row_indices, int *S_col_ptr, 
                                 int *Sk_max_sizes, int *Rk_actual_sizes, 
                                 int *rk_idx, int *submatrix, int max_Sk_size, int max_Rk_size, int n) {
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    if (col < n) {
        int start_S = S_col_ptr[col];
        int end_S = S_col_ptr[col + 1];
        
        // Initialize rk_idx and submatrix with -1
        for (int i = 0; i < Rk_actual_sizes[col]; ++i) {
            rk_idx[col * max_Sk_size + i] = -1;
        }
        for (int i = 0; i < max_Rk_size * max_Sk_size; ++i) {
            submatrix[col * max_Rk_size * max_Sk_size + i] = -1;
        }

        int rk_count = 0;
        for (int i = start_S; i < end_S; ++i) {
            int row = S_row_indices[i];
            int start_A = A_col_ptr[row];
            int end_A = A_col_ptr[row + 1];
            for (int j = start_A; j < end_A; ++j) {
                int idx = A_row_indices[j];
                bool found = false;
                for (int k = 0; k < rk_count; ++k) {
                    if (rk_idx[col * max_Sk_size + k] == idx) {
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    rk_idx[col * max_Sk_size + rk_count] = idx;
                    rk_count++;
                }
                // Add entry to the submatrix
                int submatrix_idx = col * max_Rk_size * max_Sk_size + rk_count * max_Sk_size + (i - start_S);
                submatrix[submatrix_idx] = A_values[j];
            }
        }
        Rk_actual_sizes[col] = rk_count;
    }
}
