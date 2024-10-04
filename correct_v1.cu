#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cuda_runtime.h>
#include <stdio.h>
#include <thrust/device_vector.h>
#include <thrust/extrema.h>
#include <iomanip>
#include <math.h>

#define N 100  // Adjust as needed based on your matrices

// Function to read a matrix in CSC format from a file
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

// Kernel to compute Sk_max_sizes and Rk_max_sizes
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

// Kernel to compute sk for each k
__global__ void computeSk(int *S_row_indices, int *S_col_ptr, int numCols, int *sk, int *Sk_actual_sizes, int max_size) {
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= numCols) return;  

    int startIdx = S_col_ptr[k];
    int endIdx = S_col_ptr[k + 1];
    int count = 0;
    for (int i = 0; i < max_size; i++) {
        sk[k * max_size + i] = -1;
    }
    for (int i = startIdx; i < endIdx && count < max_size; ++i) {
        sk[k * max_size + count] = S_row_indices[i]; 
        count++;
    }
    Sk_actual_sizes[k] = count;
}

// Kernel to compute rk for each k using a bitmap
__global__ void computeRk(
    int *A_row_indices, int *A_col_ptr,
    int *sk, int *Sk_actual_sizes,
    int numCols, int *rk, int *Rk_actual_sizes,
    int Rk_max_size, int Sk_max_size)
{
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= numCols) return;

    // Use bitmap for lookup to avoid excessive memory usage
    unsigned int rowBitmap[(N + 31) / 32] = {0};

    int rkCount = 0;
    int sk_size = Sk_actual_sizes[k];

    for (int i = 0; i < sk_size; i++) {
        int colIdx = sk[k * Sk_max_size + i];
        if (colIdx == -1) continue;

        int startIdx = A_col_ptr[colIdx];
        int endIdx = A_col_ptr[colIdx + 1];
        for (int j = startIdx; j < endIdx; j++) {
            int rowIdx = A_row_indices[j];
            int arrayIndex = rowIdx / 32;
            int bitPosition = rowIdx % 32;
            unsigned int mask = 1 << bitPosition;

            if ((rowBitmap[arrayIndex] & mask) == 0 && rkCount < Rk_max_size) {
                // Set the bit
                rowBitmap[arrayIndex] |= mask;
                // Add to rk
                rk[k * Rk_max_size + rkCount] = rowIdx;
                rkCount++;
            }
        }
    }
    Rk_actual_sizes[k] = rkCount;

    // Fill remaining rk entries with -1
    for (int i = rkCount; i < Rk_max_size; i++) {
        rk[k * Rk_max_size + i] = -1;
    }
}

// Device function for QR decomposition using Modified Gram-Schmidt
__device__ void qrDecomposition(double* submatrix, double* Q, double* R, int actual_rows, int actual_cols, int padded_rows, int padded_cols) {
    for (int i = 0; i < actual_cols; i++) {
        double norm = 0.0;
        for (int j = 0; j < actual_rows; j++) {
            norm += submatrix[j * padded_cols + i] * submatrix[j * padded_cols + i];
        }
        norm = sqrt(norm);
        R[i * padded_cols + i] = norm;  

        for (int j = 0; j < actual_rows; j++) {
            Q[j * padded_cols + i] = submatrix[j * padded_cols + i] / norm;
        }

        for (int j = i + 1; j < actual_cols; j++) {
            double projection = 0.0;
            for (int k = 0; k < actual_rows; k++) {
                projection += Q[k * padded_cols + i] * submatrix[k * padded_cols + j];
            }
            R[i * padded_cols + j] = projection;  

            for (int k = 0; k < actual_rows; k++) {
                submatrix[k * padded_cols + j] -= projection * Q[k * padded_cols + i];
            }
        }
    }
}

// Device function to compute Q^T * b
__device__ void computeQtb(double* Q, double* b, double* Qtb, int Q_num_rows, int Q_num_cols, int padded_cols) {
    for (int i = 0; i < Q_num_cols; i++) {
        double sum = 0.0;
        for (int j = 0; j < Q_num_rows; j++) {
            sum += Q[j * padded_cols + i] * b[j];
        }
        Qtb[i] = sum;
    }
}

// Device function for back substitution to solve R x = Q^T b
__device__ void backSubstitution(double* R, double* Qtb, double* x_local, int n, int padded_cols) {
    for (int i = n - 1; i >= 0; --i) {
        double sum = Qtb[i];
        for (int j = i + 1; j < n; ++j) {
            sum -= R[i * padded_cols + j] * x_local[j];
        }
        x_local[i] = sum / R[i * padded_cols + i];
    }
}

// Kernel to generate the submatrix for each column k and solve the system
__global__ void generateSubmatrixAndSolve(
    int *A_row_indices, int *A_col_ptr, double *A_values,
    int *sk, int *Sk_actual_sizes, int Sk_max_size,
    int *rk, int *Rk_actual_sizes, int Rk_max_size,
    int *b_row_indices, int *b_col_ptr, double *b_values,
    double *x_solutions)
{
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= N) return;

    int Sk_size = Sk_actual_sizes[k];
    int Rk_size = Rk_actual_sizes[k];

    // This shouldn't happen !
    if (Sk_size == 0 || Rk_size == 0) {
        return;
    }

    // Allocate per-thread array for row mapping
    int row_global_to_local[N]; 
    for (int i = 0; i < N; i++) {
        row_global_to_local[i] = -1;
    }

    // Build mapping from global row indices to local indices
    for (int i = 0; i < Rk_size; i++) {
        int row_global = rk[k * Rk_max_size + i];
        row_global_to_local[row_global] = i;
    }

    // Allocate local arrays
    double submatrix_local[N*N]; 
    double Q_local[N*N];
    double R_local[N*N];
    double b_local[N];
    double Qtb_local[N];
    double x_local[N];

    // Build the submatrix
    for (int col_local = 0; col_local < Sk_size; col_local++) {
        int col_global = sk[k * Sk_max_size + col_local];
        int startIdx = A_col_ptr[col_global];
        int endIdx = A_col_ptr[col_global + 1];
        for (int idx = startIdx; idx < endIdx; idx++) {
            int row_global = A_row_indices[idx];
            int row_local = row_global_to_local[row_global];
            if (row_local != -1) {
                int submatrix_index = row_local * Sk_size + col_local;
                submatrix_local[submatrix_index] = A_values[idx];
            }
        }
    }

    // Generate the local b vector


    int b_col_start = b_col_ptr[k];
    int b_col_end = b_col_ptr[k + 1];

    for (int idx = b_col_start; idx < b_col_end; idx++) {
        int b_row = b_row_indices[idx];
        int row_local = row_global_to_local[b_row];
        if (row_local != -1) {
            b_local[row_local] = b_values[idx];
        }
    }


    // Modified Gram-Schmidt
    
    // Perform QR decomposition
    qrDecomposition(submatrix_local, Q_local, R_local, Rk_size, Sk_size, Rk_size, Sk_size);

    // Compute Q^T * b
    computeQtb(Q_local, b_local, Qtb_local, Rk_size, Sk_size, Sk_size);

    // Solve R x = Q^T b
    backSubstitution(R_local, Qtb_local, x_local, Sk_size, Sk_size);

    // Store the solution x in global memory
    int x_offset = k * Sk_max_size;
    for (int i = 0; i < Sk_size; ++i) {
        x_solutions[x_offset + i] = x_local[i];
    }
}

int main() {
    std::vector<double> values_A, values_S, values_b;
    std::vector<int> rowIndices_A, rowIndices_S, rowIndices_b;
    std::vector<int> colPointers_A, colPointers_S, colPointers_b;

    // Generate filename based on the matrix size
    std::string filename1 = "matrix_1.csc";
    std::string filename2 = "matrix_2.csc"; 
    std::string filename3 = "matrix_1.csc"; 
    // Read the matrix from a file for both A and S
    readCSCFromFile(filename2, values_A, rowIndices_A, colPointers_A);
    readCSCFromFile(filename1, values_S, rowIndices_S, colPointers_S);
    readCSCFromFile(filename3, values_b, rowIndices_b, colPointers_b);

    // Allocate device memory
    int *d_rowIndices_A, *d_colPointers_A;
    int *d_rowIndices_S, *d_colPointers_S;
    int *d_rowIndices_b, *d_colPointers_b;
    int *d_Sk_max_sizes, *d_Rk_max_sizes, *d_Rk_actual_sizes, *d_Sk_actual_sizes;
    double *d_values_A;
    double *d_values_b;

    cudaMalloc(&d_rowIndices_A, rowIndices_A.size() * sizeof(int));
    cudaMalloc(&d_colPointers_A, colPointers_A.size() * sizeof(int));
    cudaMalloc(&d_rowIndices_S, rowIndices_S.size() * sizeof(int));
    cudaMalloc(&d_colPointers_S, colPointers_S.size() * sizeof(int));
    cudaMalloc(&d_rowIndices_b, rowIndices_b.size() * sizeof(int));
    cudaMalloc(&d_colPointers_b, colPointers_b.size() * sizeof(int));

    cudaMalloc(&d_Sk_max_sizes, N * sizeof(int));
    cudaMalloc(&d_Rk_max_sizes, N * sizeof(int));
    cudaMalloc(&d_Rk_actual_sizes, N * sizeof(int));  
    cudaMalloc(&d_Sk_actual_sizes, N * sizeof(int));  

    cudaMalloc(&d_values_A, values_A.size() * sizeof(double));  
    cudaMalloc(&d_values_b, values_b.size() * sizeof(double));  

    cudaMemcpy(d_rowIndices_A, rowIndices_A.data(), rowIndices_A.size() * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_colPointers_A, colPointers_A.data(), colPointers_A.size() * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_rowIndices_S, rowIndices_S.data(), rowIndices_S.size() * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_colPointers_S, colPointers_S.data(), colPointers_S.size() * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_rowIndices_b, rowIndices_b.data(), rowIndices_b.size() * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_colPointers_b, colPointers_b.data(), colPointers_b.size() * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_values_A, values_A.data(), values_A.size() * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_values_b, values_b.data(), values_b.size() * sizeof(double), cudaMemcpyHostToDevice);

    // Launch the getSkAndRkMaxSize kernel
    int blockSize = 256;
    int numBlocks = (N + blockSize - 1) / blockSize;
    getSkAndRkMaxSize<<<numBlocks, blockSize>>>(d_rowIndices_S, d_colPointers_S, d_colPointers_A, d_Sk_max_sizes, d_Rk_max_sizes, N);

    // Use Thrust to find the maximum values in d_Sk_max_sizes and d_Rk_max_sizes
    thrust::device_ptr<int> d_Sk_max_sizes_ptr(d_Sk_max_sizes);
    thrust::device_ptr<int> d_Rk_max_sizes_ptr(d_Rk_max_sizes);

    int Sk_max_size = *(thrust::max_element(d_Sk_max_sizes_ptr, d_Sk_max_sizes_ptr + N));
    int Rk_max_size = *(thrust::max_element(d_Rk_max_sizes_ptr, d_Rk_max_sizes_ptr + N));

    // Allocate memory for sk, rk_idx, submatrix, Q, and R
    int *d_sk;
    double *d_submatrix, *d_b_matrix, *d_Q, *d_R;
    cudaMalloc(&d_sk, N * Sk_max_size * sizeof(int)); 
    cudaMalloc(&d_submatrix, N * Rk_max_size * Sk_max_size * sizeof(double));
    cudaMalloc(&d_b_matrix, N * Rk_max_size  * sizeof(double));
    cudaMalloc(&d_Q, N * Rk_max_size * Sk_max_size * sizeof(double));
    cudaMalloc(&d_R, N * Sk_max_size * Sk_max_size * sizeof(double));

    // Allocate space for rk_idx, Sk_actual_sizes, Rk_actual_sizes
    int *d_rk_idx;
    double *d_x_solutions;
    double *d_Qtb;
    cudaMalloc(&d_rk_idx, N * Rk_max_size * sizeof(int));
    cudaMalloc(&d_x_solutions, N * Sk_max_size  * sizeof(double));
    cudaMalloc(&d_Qtb, N * Sk_max_size * sizeof(double));

    // Compute Sk
    computeSk<<<numBlocks, blockSize>>>(d_rowIndices_S, d_colPointers_S, N, d_sk, d_Sk_actual_sizes, Sk_max_size);

    // Compute Rk
    computeRk<<<numBlocks, blockSize>>>(
        d_rowIndices_A, d_colPointers_A,
        d_sk, d_Sk_actual_sizes,
        N, d_rk_idx, d_Rk_actual_sizes,
        Rk_max_size, Sk_max_size);

    // Generate Submatrices and solve the systems
    generateSubmatrixAndSolve<<<numBlocks, blockSize>>>(
        d_rowIndices_A, d_colPointers_A, d_values_A,
        d_sk, d_Sk_actual_sizes, Sk_max_size,
        d_rk_idx, d_Rk_actual_sizes, Rk_max_size,
        d_rowIndices_b, d_colPointers_b, d_values_b,
        d_x_solutions);

    // Copy submatrices from device to host
    std::vector<double> h_submatrix(N * Rk_max_size * Sk_max_size);
    cudaMemcpy(h_submatrix.data(), d_submatrix, N * Rk_max_size * Sk_max_size * sizeof(double), cudaMemcpyDeviceToHost);

    // Copy Sk_actual_sizes and Rk_actual_sizes to host
    std::vector<int> h_Sk_actual_sizes(N);
    std::vector<int> h_Rk_actual_sizes(N);
    cudaMemcpy(h_Sk_actual_sizes.data(), d_Sk_actual_sizes, N * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_Rk_actual_sizes.data(), d_Rk_actual_sizes, N * sizeof(int), cudaMemcpyDeviceToHost);

    // Copy sk and rk_idx to host
    std::vector<int> h_sk(N * Sk_max_size);
    std::vector<int> h_rk_idx(N * Rk_max_size);
    cudaMemcpy(h_sk.data(), d_sk, N * Sk_max_size * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_rk_idx.data(), d_rk_idx, N * Rk_max_size * sizeof(int), cudaMemcpyDeviceToHost);

    // Copy x_solutions from device to host
    std::vector<double> h_x_solutions(N * Sk_max_size);
    cudaMemcpy(h_x_solutions.data(), d_x_solutions, N * Sk_max_size * sizeof(double), cudaMemcpyDeviceToHost);

    // Print the solutions
    for (int k = 0; k < N; ++k) {
        int Sk_size = h_Sk_actual_sizes[k];
        if (Sk_size == 0) {
            continue;
        }
        std::cout << "Solution x for column " << k << ":" << std::endl;
        for (int i = 0; i < Sk_size; ++i) {
            int col_idx = h_sk[k * Sk_max_size + i];
            double value = h_x_solutions[k * Sk_max_size + i];
            std::cout << "x[" << col_idx << "] = " << value << std::endl;
        }
        std::cout << "-----------------------------" << std::endl;
    }

    // Free device memory
    cudaFree(d_rowIndices_A);
    cudaFree(d_colPointers_A);
    cudaFree(d_rowIndices_S);
    cudaFree(d_colPointers_S);
    cudaFree(d_rowIndices_b);
    cudaFree(d_colPointers_b);
    cudaFree(d_values_b);
    cudaFree(d_Sk_max_sizes);
    cudaFree(d_Rk_max_sizes);
    cudaFree(d_Rk_actual_sizes);
    cudaFree(d_Sk_actual_sizes);
    cudaFree(d_sk);
    cudaFree(d_rk_idx);
    cudaFree(d_values_A);
    cudaFree(d_x_solutions);

    return 0;
}
