__global__ void getSkAndRkMaxSize(
    int *S_row_indices, int *S_col_ptr, int *A_col_ptr, 
    int *Sk_max_sizes, int *Rk_max_sizes, int n, 
    int *Sk_max_global, int *Rk_max_global) {
    
    extern __shared__ int shared_mem[]; // Dynamically allocated shared memory
    int *Sk_max_shared = shared_mem;    // Shared memory for Sk_max_size
    int *Rk_max_shared = &shared_mem[blockDim.x]; // Shared memory for Rk_max_size

    int col = blockIdx.x * blockDim.x + threadIdx.x;
    
    // Initialize shared memory for reduction
    Sk_max_shared[threadIdx.x] = 0;
    Rk_max_shared[threadIdx.x] = 0;
    __syncthreads();

    if (col < n) {
        int Sk_max_size = 0;
        int Rk_max_size = 0;

        int start_S = S_col_ptr[col];
        int end_S = S_col_ptr[col + 1];
        Sk_max_size = end_S - start_S;  // Number of non-zero entries in column col of S

        for (int i = start_S; i < end_S; ++i) {
            int row = S_row_indices[i];
            int start_A = A_col_ptr[row];
            int end_A = A_col_ptr[row + 1];
            Rk_max_size += (end_A - start_A);
        }

        Sk_max_sizes[col] = Sk_max_size;
        Rk_max_sizes[col] = Rk_max_size;

        // Update shared memory
        Sk_max_shared[threadIdx.x] = Sk_max_size;
        Rk_max_shared[threadIdx.x] = Rk_max_size;
    }

    __syncthreads();

    // Perform parallel reduction within the block to find the max in shared memory
    for (int stride = blockDim.x / 2; stride > 0; stride >>= 1) {
        if (threadIdx.x < stride) {
            Sk_max_shared[threadIdx.x] = max(Sk_max_shared[threadIdx.x], Sk_max_shared[threadIdx.x + stride]);
            Rk_max_shared[threadIdx.x] = max(Rk_max_shared[threadIdx.x], Rk_max_shared[threadIdx.x + stride]);
        }
        __syncthreads();
    }

    // Atomic update to global maximums
    if (threadIdx.x == 0) {
        atomicMax(Sk_max_global, Sk_max_shared[0]);
        atomicMax(Rk_max_global, Rk_max_shared[0]);
    }
}


int blockSize = 256;
int numBlocks = (N + blockSize - 1) / blockSize;
size_t sharedMemSize = 2 * blockSize * sizeof(int); // Shared memory for Sk and Rk maxes

int h_Sk_max_global = 0, h_Rk_max_global = 0;
int *d_Sk_max_global, *d_Rk_max_global;
cudaMalloc(&d_Sk_max_global, sizeof(int));
cudaMalloc(&d_Rk_max_global, sizeof(int));
cudaMemcpy(d_Sk_max_global, &h_Sk_max_global, sizeof(int), cudaMemcpyHostToDevice);
cudaMemcpy(d_Rk_max_global, &h_Rk_max_global, sizeof(int), cudaMemcpyHostToDevice);

getSkAndRkMaxSize<<<numBlocks, blockSize, sharedMemSize>>>(
    d_rowIndices_S, d_colPointers_S, d_colPointers_A, 
    d_Sk_max_sizes, d_Rk_max_sizes, N, 
    d_Sk_max_global, d_Rk_max_global);

cudaMemcpy(&h_Sk_max_global, d_Sk_max_global, sizeof(int), cudaMemcpyDeviceToHost);
cudaMemcpy(&h_Rk_max_global, d_Rk_max_global, sizeof(int), cudaMemcpyDeviceToHost);

// Now h_Sk_max_global and h_Rk_max_global hold the max values
