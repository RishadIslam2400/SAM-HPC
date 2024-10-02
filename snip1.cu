#include <cuda_runtime.h>

/// provided 2 ways to get the max size of Sk and Rk, need to evaluate which is better, this is just for avoid out of kernel reduction, other optimization can be applied orthognally

// version 1: on stack reduction
__global__ void getSkAndRkMaxSize(
    int *S_row_indices, int *S_col_ptr, int *A_col_ptr, 
    int *Sk_max_sizes, int *Rk_max_sizes, int n, 
    int *Sk_max_global, int *Rk_max_global) {

    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int Sk_max_size = 0;
    int Rk_max_size = 0;

    if (col < n) {
        int start_S = S_col_ptr[col];
        int end_S = S_col_ptr[col + 1];
        Sk_max_size = end_S - start_S;  // Number of non-zero entries in column col of S

        for (int i = start_S; i < end_S; ++i) {
            int row = S_row_indices[i];
            int start_A = A_col_ptr[row];
            int end_A = A_col_ptr[row + 1];
            Rk_max_size += (end_A - start_A);
        }
    }

    // Synchronize to ensure all threads have computed their values
    __syncthreads();

    // Use atomicMax to update the global maximum directly from the on-stack values
    if (col < n) {
        atomicMax(Sk_max_global, Sk_max_size);
        atomicMax(Rk_max_global, Rk_max_size);
    }
}

// version 2: shared memory reduction
__global__ void getSkAndRkMaxSize(
    int *S_row_indices, int *S_col_ptr, int *A_col_ptr, 
    int n, int *Sk_max_global, int *Rk_max_global) {

    extern __shared__ int shared_mem[]; // Dynamically allocated shared memory
    int *Sk_max_shared = shared_mem;    // Shared memory for Sk_max_size
    int *Rk_max_shared = &shared_mem[blockDim.x]; // Shared memory for Rk_max_size

    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int Sk_max_size = 0;
    int Rk_max_size = 0;

    if (col < n) {
        int start_S = S_col_ptr[col];
        int end_S = S_col_ptr[col + 1];
        Sk_max_size = end_S - start_S;  // Number of non-zero entries in column col of S

        for (int i = start_S; i < end_S; ++i) {
            int row = S_row_indices[i];
            int start_A = A_col_ptr[row];
            int end_A = A_col_ptr[row + 1];
            Rk_max_size += (end_A - start_A);
        }
    }

    // Store the results in shared memory
    Sk_max_shared[threadIdx.x] = Sk_max_size;
    Rk_max_shared[threadIdx.x] = Rk_max_size;

    // Synchronize to make sure all threads have written to shared memory
    __syncthreads();

    // Perform block-level reduction in shared memory
    for (int stride = blockDim.x / 2; stride > 0; stride >>= 1) {
        if (threadIdx.x < stride) {
            Sk_max_shared[threadIdx.x] = max(Sk_max_shared[threadIdx.x], Sk_max_shared[threadIdx.x + stride]);
            Rk_max_shared[threadIdx.x] = max(Rk_max_shared[threadIdx.x], Rk_max_shared[threadIdx.x + stride]);
        }
        __syncthreads(); // Ensure all threads complete the reduction step
    }

    // The first thread in each block stores the block-level max to global memory
    if (threadIdx.x == 0) {
        Sk_max_global[blockIdx.x] = Sk_max_shared[0];
        Rk_max_global[blockIdx.x] = Rk_max_shared[0];
    }
}
