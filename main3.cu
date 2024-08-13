#include <cuda_runtime.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdio.h>
#include <string>
#include <thrust/device_vector.h>
#include <thrust/extrema.h>
#include <vector>

// Device function for binary search
__device__ int binarySearch(int *arr, int low, int high, int target) {
  while (low <= high) {
    int mid = low + (high - low) / 2;

    // Check if the target is present at mid
    if (arr[mid] == target) {
      return mid;
    }

    // If target is greater, ignore the left half
    if (arr[mid] < target) {
      low = mid + 1;
    }
    // If target is smaller, ignore the right half
    else {
      high = mid - 1;
    }
  }
  // If the target is not present in the array
  return -1;
}

// Device function to swap elements
__device__ void swap(int *arr, int i, int j) {
  int temp = arr[i];
  arr[i] = arr[j];
  arr[j] = temp;
}

// Device function to partition the array
__device__ int partition(int *arr, int low, int high) {
  int pivot = arr[high];
  int i = low - 1;
  for (int j = low; j < high; ++j) {
    if (arr[j] < pivot) {
      ++i;
      swap(arr, i, j);
    }
  }
  swap(arr, i + 1, high);
  return i + 1;
}

// Device function for quicksort
__device__ void quicksort(int *arr, int low, int high) {
  if (low < high) {
    int pi = partition(arr, low, high);
    quicksort(arr, low, pi - 1);
    quicksort(arr, pi + 1, high);
  }
}

void readCSCFromFile(
    const std::string &filename,
    std::vector<double> &values, // each value is the non-zero value
    std::vector<int>
        &rowIndices, // each value is row index of the non-zero value
    std::vector<int>
        &colPointers) { // each value index into rowIndices vector, point to the
                        // first non-zero value in that column
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

void printCSCArrays(const std::vector<double> &values,
                    const std::vector<int> &rowIndices,
                    const std::vector<int> &colPointers) {
  std::cout << "Values Array: ";
  for (const auto &val : values) {
    std::cout << val << " ";
  }
  std::cout << std::endl;

  std::cout << "Row Indices Array: ";
  for (const auto &row : rowIndices) {
    std::cout << row << " ";
  }
  std::cout << std::endl;

  std::cout << "Column Pointers Array: ";
  for (const auto &col : colPointers) {
    std::cout << col << " ";
  }
  std::cout << std::endl;
}

__global__ void getSkAndRkMaxSize(int *S_row_indices, int *S_col_ptr,
                                  int *A_col_ptr, int *Sk_max_sizes,
                                  int *Rk_max_sizes, int n) {
  // iterate over each column of S
  int col = blockIdx.x * blockDim.x + threadIdx.x;
  if (col < n) {
    int Sk_max_size = 0;
    int Rk_max_size = 0;
    int start_S = S_col_ptr[col];
    int end_S = S_col_ptr[col + 1];
    Sk_max_size =
        end_S - start_S; // Number of non-zero entries in column `col` of S

    for (int i = start_S; i < end_S; ++i) {
      int row = S_row_indices[i]; // non-zero entry's row in column `col` of S;
                                  // row also index into A's column pointers;
      int start_A = A_col_ptr[row]; // start_A index into A's row indices
      int end_A = A_col_ptr[row + 1];
      Rk_max_size += (end_A - start_A);
    }

    Sk_max_sizes[col] = Sk_max_size;
    Rk_max_sizes[col] = Rk_max_size;
  }
}

__global__ void extractSubMatrix(double *A_values, int *A_row_indices,
                                 int *A_col_ptr, int *S_row_indices,
                                 int *S_col_ptr, int *Sk_max_sizes,
                                 int *Rk_actual_sizes, int *rk_idx,
                                 double *submatrix, int max_Sk_size,
                                 int max_Rk_size, int n) {
  int col = blockIdx.x * blockDim.x + threadIdx.x;
  if (col < n) {
    int start_S = S_col_ptr[col];
    int end_S = S_col_ptr[col + 1];
    int Sk_size = end_S - start_S; // Number of non-zero entries in column `col`
                                   // of S; also index into A's column pointers

    // Initialize rk_idx with -1
    for (int i = 0; i < max_Rk_size; ++i) {
      rk_idx[col * max_Rk_size + i] = -1;
    }
    for (int i = 0; i < max_Rk_size * max_Sk_size; ++i) {
      submatrix[col * max_Rk_size * max_Sk_size + i] =
          0.0; // Initialize submatrix with zeros
    }

    int rk_count = 0; // number of rows of A that contains non-zero entries
    for (int i = start_S; i < end_S; ++i) { // iterate over col of A
      int row = S_row_indices[i];   // Row index in column `col` of S also index
                                    // into A's col pointers
      int start_A = A_col_ptr[row]; // start_A index into A's row indices:
      int end_A = A_col_ptr[row + 1];
      for (int j = start_A; j < end_A;
           ++j) { // process the column[row] of A, A[start_A][row] ...
                  // A[end_A][row]
        int idx = A_row_indices[j];
        bool found = false;
        for (int k = 0; k < rk_count; ++k) {
          if (rk_idx[col * max_Rk_size + k] == idx) {
            found = true;
            break;
          }
        }
        if (!found) {
          // update rk_idx to include idx and increment rk_count
          rk_idx[col * max_Rk_size + rk_count] = idx;
          rk_count++;
        }
      }
    }
    Rk_actual_sizes[col] = rk_count;
    // Sort rk_idx
    quicksort(rk_idx, col * max_Rk_size, col * max_Rk_size + rk_count - 1);
    // generate plain submatrix
    for (int i = start_S; i < end_S; ++i) {
      int row = S_row_indices[i];
      int start_A = A_col_ptr[row];
      int end_A = A_col_ptr[row + 1];
      for (int j = start_A; j < end_A; ++j) {
        int idx = A_row_indices[j];
        int submatrix_idx =
            col * max_Rk_size * max_Sk_size +
            binarySearch(rk_idx, col * max_Rk_size,
                         col * max_Rk_size + rk_count - 1, idx) *
                max_Sk_size +
            (i - start_S);
        submatrix[submatrix_idx] = A_values[j];
      }
    }
  }
}

int main() {
  int N = 10; // Example size, adjust as needed

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
  int *d_Sk_max_sizes, *d_Rk_max_sizes, *d_Rk_actual_sizes;
  double *d_values_A;

  cudaMalloc(&d_rowIndices_A, rowIndices_A.size() * sizeof(int));
  cudaMalloc(&d_colPointers_A, colPointers_A.size() * sizeof(int));
  cudaMalloc(&d_rowIndices_S, rowIndices_S.size() * sizeof(int));
  cudaMalloc(&d_colPointers_S, colPointers_S.size() * sizeof(int));
  cudaMalloc(&d_Sk_max_sizes, N * sizeof(int));
  cudaMalloc(&d_Rk_max_sizes, N * sizeof(int));
  cudaMalloc(&d_Rk_actual_sizes,
             N * sizeof(int)); // Allocate memory for actual sizes
  cudaMalloc(&d_values_A,
             values_A.size() * sizeof(double)); // Allocate memory for A values

  cudaMemcpy(d_rowIndices_A, rowIndices_A.data(),
             rowIndices_A.size() * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_colPointers_A, colPointers_A.data(),
             colPointers_A.size() * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_rowIndices_S, rowIndices_S.data(),
             rowIndices_S.size() * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_colPointers_S, colPointers_S.data(),
             colPointers_S.size() * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_values_A, values_A.data(), values_A.size() * sizeof(double),
             cudaMemcpyHostToDevice);

  // Launch the getSkAndRkMaxSize kernel
  int blockSize = 256;
  int numBlocks = (N + blockSize - 1) / blockSize;
  getSkAndRkMaxSize<<<numBlocks, blockSize>>>(d_rowIndices_S, d_colPointers_S,
                                              d_colPointers_A, d_Sk_max_sizes,
                                              d_Rk_max_sizes, N);

  // Use Thrust to find the maximum values in d_Sk_max_sizes and d_Rk_max_sizes
  thrust::device_ptr<int> d_Sk_max_sizes_ptr(d_Sk_max_sizes);
  thrust::device_ptr<int> d_Rk_max_sizes_ptr(d_Rk_max_sizes);

  int Sk_max_size =
      *(thrust::max_element(d_Sk_max_sizes_ptr, d_Sk_max_sizes_ptr + N));
  int Rk_max_size =
      *(thrust::max_element(d_Rk_max_sizes_ptr, d_Rk_max_sizes_ptr + N));

  // Allocate memory for rk_idx and submatrix
  int *d_rk_idx;
  double *d_submatrix;
  cudaMalloc(&d_rk_idx, N * Rk_max_size * sizeof(int));
  cudaMalloc(&d_submatrix, N * Rk_max_size * Sk_max_size * sizeof(double));

  // Launch the extractSubMatrix kernel
  extractSubMatrix<<<numBlocks, blockSize>>>(
      d_values_A, d_rowIndices_A, d_colPointers_A, d_rowIndices_S,
      d_colPointers_S, d_Sk_max_sizes, d_Rk_actual_sizes, d_rk_idx, d_submatrix,
      Sk_max_size, Rk_max_size, N);

  // Copy rk_idx, d_Rk_actual_sizes, and submatrix back to host
  std::vector<int> h_rk_idx(N * Rk_max_size);
  std::vector<int> h_Rk_actual_sizes(N);
  std::vector<double> h_submatrix(N * Rk_max_size * Sk_max_size);

  cudaMemcpy(h_rk_idx.data(), d_rk_idx, N * Rk_max_size * sizeof(int),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(h_Rk_actual_sizes.data(), d_Rk_actual_sizes, N * sizeof(int),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(h_submatrix.data(), d_submatrix,
             N * Rk_max_size * Sk_max_size * sizeof(double),
             cudaMemcpyDeviceToHost);

  // Print the Rk_actual_sizes results
  std::cout << "Rk_actual_sizes:" << std::endl;
  for (int i = 0; i < N; ++i) {
    std::cout << "Column " << i << ": " << h_Rk_actual_sizes[i] << std::endl;
  }

  // Print the rk_idx results
  for (int i = 0; i < N; ++i) {
    std::cout << "rk_idx for column " << i << ": ";
    for (int j = 0; j < Rk_max_size; ++j) {
      if (h_rk_idx[i * Rk_max_size + j] == -1)
        break; // End of valid entries
      std::cout << h_rk_idx[i * Rk_max_size + j] << " ";
    }
    std::cout << std::endl;
  }

  // Print the submatrix results
  std::cout << "Submatrix:" << std::endl;
  for (int i = 0; i < N; ++i) {
    std::cout << "Submatrix for column " << i << ": " << std::endl;
    for (int j = 0; j < Rk_max_size; ++j) {
      for (int k = 0; k < Sk_max_size; ++k) {
        double value =
            h_submatrix[i * Rk_max_size * Sk_max_size + j * Sk_max_size + k];
        std::cout << std::setw(5) << value << " ";
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
  cudaFree(d_values_A);

  return 0;
}