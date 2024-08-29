#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cuda_runtime.h>
#include <stdio.h>
#include <thrust/device_vector.h>
#include <thrust/extrema.h>
#include <iomanip>

#define N 10



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

std::vector<double> incrementValues(const std::vector<double>& values) {
    std::vector<double> incrementedValues(values.size());
    for (size_t i = 0; i < values.size(); ++i) {
        incrementedValues[i] = values[i] + 10.0;
    }
    return incrementedValues;
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


__global__ void computeRk(int *A_row_indices, int *A_col_ptr, int *sk, int numCols, int *rk, int max_size) {
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= numCols) return;

    int lookup[N]={0};  
    for (int i = 0; i < max_size; i++) {
        rk[k * max_size + i] = -1;
    }

    int rkCount = 0;  
    for (int i = 0; i < max_size; i++) {
        int colIdx = sk[k * max_size + i];
        if (colIdx == -1) continue;  

        int startIdx = A_col_ptr[colIdx];
        int endIdx = A_col_ptr[colIdx + 1];
        for (int j = startIdx; j < endIdx; j++) {
            int rowIdx = A_row_indices[j];
            if (lookup[rowIdx] == 0 && rkCount < max_size) {  
                rk[k * max_size + rkCount] = rowIdx;
                rkCount++;
                lookup[rowIdx] = 1;  
            }
        }
    }

}


__device__ double vectorNorm(double* v, int length) {
    double norm = 0.0;
    for (int i = 0; i < length; i++) {
        norm += v[i] * v[i];
    }
    return sqrt(norm);
}


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

__device__ void computeQtb(double* Q, double* b, double* Qtb, int Q_num_rows, int Q_num_cols, int max_Sk_size) {
    
    int Qtb_num_rows = Q_num_cols;
    for (int i = 0; i < Q_num_cols; i++) {
        double sum = 0.0;
        for (int j = 0; j < Q_num_rows; j++) {
            sum += Q[j * max_Sk_size + i] * b[j];
        }
        Qtb[i] = sum;
    }
}

__device__ void computeQtb1(double* Q, double* b, double* Qtb, int Q_num_rows, int Q_num_cols, int max_Sk_size) {
    
    int Qtb_num_rows = Q_num_cols;
    for (int i = 0; i < Q_num_cols; i++) {
        double sum = 0.0;
        for (int j = 0; j < Q_num_rows; j++) {
            // printf("Q[%d][%d] = %f, b[%d] = %f\n", j, i, Q[j * max_Sk_size + i], j, b[j]);
            sum += Q[j * max_Sk_size + i] * b[j];
        }
        Qtb[i] = sum;
    }
}



__device__ void backSubstitution(double* R, double* Qtb, double* x_matrix, int n, int col, int max_Sk_size) {
    int x_start = col * max_Sk_size;
    for (int i = 0; i < n; i++) {
        x_matrix[x_start + i] = Qtb[i]/ R[i * max_Sk_size + i];
    }

    // int x_index = col * max_Rk_size; // Calculate the starting index for the solution in x_matrix
    // for (int i = n - 1; i >= 0; i--) {
    //     double sum = 0.0;
    //     for (int j = i + 1; j < n; j++) {
    //         sum += R[i * n + j] * x_matrix[x_index + j]; // Use the correctly offset x values
    //     }
    //     x_matrix[x_index + i] = (Qtb[i] - sum) / R[i * n + i];  // Compute and store directly in x_matrix
    // }
}




// Kernel that performs submatrix extraction followed by QR decomposition
__global__ void kernel1(double *A_values, double *b_values,
                                 int *A_row_indices, int *A_col_ptr, 
                                 int *S_row_indices, int *S_col_ptr, 
                                 int *b_row_indices, int *b_col_ptr, 
                                 int *Rk_actual_sizes, int *Sk_actual_sizes,
                                 int *rk_idx, double *submatrix, double *b_matrix,
                                 double *Q, double *R, 
                                 double *Qtb, double *x_matrix,
                                 int max_Sk_size, int max_Rk_size, int n) {
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    if (col < n) {
        int start_S = S_col_ptr[col];
        int end_S = S_col_ptr[col + 1];

        // Initialize rk_idx and submatrix with -1
        for (int i = 0; i < max_Rk_size; ++i) {
            rk_idx[col * max_Rk_size + i] = -1;
        }
        for (int i = 0; i < max_Rk_size * max_Sk_size; ++i) {
            submatrix[col * max_Rk_size * max_Sk_size + i] = 0.0;  // Initialize submatrix with zeros
        }
        for (int i = 0; i < max_Rk_size; ++i) {
            b_matrix[col * max_Rk_size + i] = 0.0;
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
                    if (rk_idx[col * max_Rk_size + k] == idx) {
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    rk_idx[col * max_Rk_size + rk_count] = idx;
                    // Add entry to the submatrix
                    for (int l = start_S; l < end_S; ++l) {
                        int sub_row = S_row_indices[l];
                        int sub_start_A = A_col_ptr[sub_row];
                        int sub_end_A = A_col_ptr[sub_row + 1];
                        bool value_found = false;
                        for (int m = sub_start_A; m < sub_end_A; ++m) {
                            if (A_row_indices[m] == idx) {
                                int submatrix_idx = col * max_Rk_size * max_Sk_size + rk_count * max_Sk_size + (l - start_S);
                                submatrix[submatrix_idx] = A_values[m];
                                value_found = true;
                                break;
                            }
                        }
                        if (!value_found) {
                            int submatrix_idx = col * max_Rk_size * max_Sk_size + rk_count * max_Sk_size + (l - start_S);
                            submatrix[submatrix_idx] = 0.0;
                        }
                    }
                    
                    int start_b = b_col_ptr[col];
                    int end_b = b_col_ptr[col + 1];
                    bool row_idx_found = false;
                    for (int l = start_b; l < end_b; ++l) {
                        if (b_row_indices[l] == idx) {
                            b_matrix[col * max_Rk_size + rk_count] = b_values[l];
                            row_idx_found = true;
                            break;
                        }
                    }
                    if (!row_idx_found) {
                        b_matrix[col * max_Rk_size + rk_count] = 0.0;
                    }
                    rk_count++;
                }
            }
        }
        Rk_actual_sizes[col] = rk_count;

        qrDecomposition(&submatrix[col * max_Rk_size * max_Sk_size],
                        &Q[col * max_Rk_size * max_Sk_size],
                        &R[col * max_Sk_size * max_Sk_size],
                        rk_count, end_S - start_S, max_Rk_size, max_Sk_size);


        if (col==0){
            computeQtb1(&Q[col * max_Rk_size * max_Sk_size], &b_matrix[col * max_Rk_size],
                   &Qtb[col * max_Sk_size], Rk_actual_sizes[col], Sk_actual_sizes[col], max_Sk_size);
        }else{
            computeQtb(&Q[col * max_Rk_size * max_Sk_size], &b_matrix[col * max_Rk_size],
                   &Qtb[col * max_Sk_size], Rk_actual_sizes[col], Sk_actual_sizes[col], max_Sk_size);
        }


        backSubstitution(&R[col * max_Sk_size * max_Sk_size], &Qtb[col * max_Sk_size], x_matrix, Sk_actual_sizes[col], col, max_Sk_size);

        if (col == 0)
        {

            printf("Rk_actual_sizes[%d] = %d, Sk_actual_sizes[%d] = %d\n", col, Rk_actual_sizes[col], col, Sk_actual_sizes[col]);
            // Print Q matrix for column 0
            printf("Matrix Q (column 0):\n");
            for (int i = 0; i < max_Rk_size; i++) {
                for (int j = 0; j < max_Sk_size; j++) {
                    printf("%f ", Q[col * max_Rk_size * max_Sk_size + i * max_Sk_size + j]);
                }
                printf("\n");
            }
            
            // Print b vector for column 0
            printf("Vector b (column 0):\n");
            for (int i = 0; i < max_Rk_size; i++) {
                printf("%f\n", b_matrix[col * max_Sk_size + i]);
            }
            // Print R matrix for column 1
            printf("Matrix R (column 0):\n");
            for (int i = 0; i < max_Sk_size; i++) {
                for (int j = 0; j < max_Sk_size; j++) {
                    printf("%f ", R[col * max_Sk_size * max_Sk_size + i * max_Sk_size + j]);
                }
                printf("\n");
            }
            
            // Print Qtb vector for column 1
            printf("Vector Qtb (column 0):\n");
            for (int i = 0; i < Sk_actual_sizes[col]; i++) {
                printf("%f\n", Qtb[col * max_Sk_size + i]);
            }
        }
    }
}







int main() {
    std::vector<double> values_A, values_S, values_b;
    std::vector<int> rowIndices_A, rowIndices_S, rowIndices_b;
    std::vector<int> colPointers_A, colPointers_S, colPointers_b;

    // Generate filename based on the matrix size
    std::string filename1 = "matrix_" + std::to_string(N) + ".csc";
    std::string filename2 = "matrixA_" + std::to_string(N) + ".csc"; 
    std::string filename3 = "matrixb_" + std::to_string(N) + ".csc"; 
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
    cudaMalloc(&d_rk_idx, N * Rk_max_size * sizeof(int));

    double *d_x_solutions;
    cudaMalloc(&d_x_solutions, N * Sk_max_size  * sizeof(double));

    double *d_Qtb;
    cudaMalloc(&d_Qtb, N * Sk_max_size * sizeof(double));

    computeSk<<<numBlocks, blockSize>>>(d_rowIndices_S, d_colPointers_S, N, d_sk, d_Sk_actual_sizes, Sk_max_size);

    kernel1<<<numBlocks, blockSize>>>(d_values_A, d_values_b,
                                               d_rowIndices_A, d_colPointers_A,
                                               d_rowIndices_S, d_colPointers_S,
                                               d_rowIndices_b, d_colPointers_b,
                                               d_Rk_actual_sizes, d_Sk_actual_sizes,
                                               d_rk_idx, d_submatrix, d_b_matrix,
                                               d_Q, d_R,
                                               d_Qtb, d_x_solutions,
                                               Sk_max_size, Rk_max_size, N);

    // Copy Q and R matrices from device to host
    std::vector<double> h_Q(N * Rk_max_size * Sk_max_size);
    std::vector<double> h_R(N * Rk_max_size * Sk_max_size);
    std::vector<double> h_submatrix(N * Rk_max_size * Sk_max_size);
    std::vector<double> h_x_solutions(N * Sk_max_size);
    std::vector<double> h_Qtb(N * Sk_max_size);


    cudaMemcpy(h_Q.data(), d_Q, N * Rk_max_size * Sk_max_size * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_R.data(), d_R, N * Sk_max_size * Sk_max_size * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_submatrix.data(), d_submatrix, N * Rk_max_size * Sk_max_size * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_x_solutions.data(), d_x_solutions, N * Sk_max_size * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_Qtb.data(), d_Qtb, N * Sk_max_size * sizeof(double), cudaMemcpyDeviceToHost);


    for (int i = 0; i < N; ++i) {
        std::cout << "Q matrix for column " << i << ":\n";
        for (int row = 0; row < Rk_max_size; ++row) {
            for (int col = 0; col < Sk_max_size; ++col) {
                std::cout << std::setw(5) << h_Q[i * Rk_max_size * Sk_max_size + row * Sk_max_size + col] << " ";
            }
            std::cout << std::endl;
        }

        std::cout << "R matrix for column " << i << ":\n";
        for (int row = 0; row < Sk_max_size; ++row) {
            for (int col = 0; col < Sk_max_size; ++col) {
                std::cout << std::setw(5) << h_R[i * Sk_max_size * Sk_max_size + row * Sk_max_size + col] << " ";
            }
            std::cout << std::endl;
        }
    }

    // Print the submatrix results
    std::cout << "Submatrix:" << std::endl;
    for (int i = 0; i < N; ++i) {
        std::cout << "Submatrix for column " << i << ": " << std::endl;
        for (int j = 0; j < Rk_max_size; ++j) {
            for (int k = 0; k < Sk_max_size; ++k) {
                double value = h_submatrix[i * Rk_max_size * Sk_max_size + j * Sk_max_size + k];
                std::cout << std::setw(5) << value << " ";
            }
            std::cout << std::endl;
        }
    }

    for (int i = 0; i < N; ++i) {
        std::cout << "Qtb " << i << ": " << std::endl;
        for (int j = 0; j < Sk_max_size; ++j) {
            double value = h_Qtb[i * Sk_max_size + j];
            std::cout << std::setw(5) << value << " ";
        }
        std::cout << std::endl;
    }



    for (int i = 0; i < N; ++i) {
        std::cout << "Solution for column " << i << ": " << std::endl;
        for (int j = 0; j < Sk_max_size; ++j) {
            double value = h_x_solutions[i * Sk_max_size + j];
            std::cout << std::setw(5) << value << " ";
        }
        std::cout << std::endl;
    }



    // Free device memory
    cudaFree(d_rowIndices_A);
    cudaFree(d_colPointers_A);
    cudaFree(d_rowIndices_S);
    cudaFree(d_colPointers_S);
    cudaFree(d_rowIndices_b);
    cudaFree(d_colPointers_b);
    cudaFree(d_Sk_max_sizes);
    cudaFree(d_Rk_max_sizes);
    cudaFree(d_Rk_actual_sizes);
    cudaFree(d_rk_idx);
    cudaFree(d_submatrix);
    cudaFree(d_values_A);
    cudaFree(d_values_b);
    cudaFree(d_Q);
    cudaFree(d_R);

    return 0;
}
