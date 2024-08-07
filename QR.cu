/* 
 *  QR Factorization Dense Linear Solver
 *
 *  filename: test_cusolver_cuda6d5.cpp
 *
 *  module add cudatoolkit
 *  compile:  nvcc -o test_cusolver_cuda6d5 test_cusolver_cuda6d5.cpp -lcublas -lcusolver
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <cuda_runtime.h>

#include <cublas_v2.h>
#include <cusolverDn.h>

void printMatrix(int m, int n, const double*A, int lda, const char* name)
{
    for(int row = 0 ; row < m ; row++){
        for(int col = 0 ; col < n ; col++){
            double Areg = A[row + col*lda];
            printf("%s(%d,%d) = %f\n", name, row+1, col+1, Areg);
        }
    } 
}

int main(int argc, char*argv[])
{
    cusolverDnHandle_t cudenseH = NULL;
    cublasHandle_t cublasH = NULL;
    cublasStatus_t cublas_status = CUBLAS_STATUS_SUCCESS;
    cusolverStatus_t cusolver_status = CUSOLVER_STATUS_SUCCESS;
    cudaError_t cudaStat1 = cudaSuccess;
    cudaError_t cudaStat2 = cudaSuccess;
    cudaError_t cudaStat3 = cudaSuccess;
    cudaError_t cudaStat4 = cudaSuccess;
    const int m = 5;
    const int n = 3;
    const int lda = m;
    const int ldb = m;
    const int nrhs = 1; // number of right hand side vectors
/*       | 1 0 2 |
 *   A = | 0 3 0 | 
 *       | 4 0 0 |
 *       | 0 0 5 |
 *       | 6 0 0 |
 *   x* = ( 1 1 2 )'
 *   b = ( 5 3 4 10 6 )'
 */
    double A[lda*n] = { 
        1.0, 0.0, 4.0, 0.0, 6.0, 
        0.0, 3.0, 0.0, 0.0, 0.0,
        2.0, 0.0, 0.0, 5.0, 0.0
        };
//    double X[ldb*nrhs] = { 1.0, 1.0, 1.0}; // exact solution
    double B[ldb*nrhs] = { 5.0, 3.0, 4.0, 10.0, 6.0};
    double XC[ldb*nrhs]; // solution matrix from GPU

    double *d_A = NULL; // linear memory of GPU
    double *d_tau = NULL; // linear memory of GPU
    double *d_B  = NULL;
    int *devInfo = NULL; // info in gpu (device copy)
    double *d_work = NULL;
    int  lwork = 0;

    int info_gpu = 0;

    const double one = 1;
    printf("A = (matlab base-1)\n");
    printMatrix(m, n, A, lda, "A");
    printf("=====\n");
    printf("B = (matlab base-1)\n");
    printMatrix(m, nrhs, B, ldb, "B");
    printf("=====\n");

// step 1: create cudense/cublas handle
    cusolver_status = cusolverDnCreate(&cudenseH);
    assert(CUSOLVER_STATUS_SUCCESS == cusolver_status);

    cublas_status = cublasCreate(&cublasH);
    assert(CUBLAS_STATUS_SUCCESS == cublas_status);

// step 2: copy A and B to device
    cudaStat1 = cudaMalloc ((void**)&d_A  , sizeof(double) * lda * n);
    cudaStat2 = cudaMalloc ((void**)&d_tau, sizeof(double) * n);
    cudaStat3 = cudaMalloc ((void**)&d_B  , sizeof(double) * ldb * nrhs);
    cudaStat4 = cudaMalloc ((void**)&devInfo, sizeof(int));
    assert(cudaSuccess == cudaStat1);
    assert(cudaSuccess == cudaStat2);
    assert(cudaSuccess == cudaStat3);
    assert(cudaSuccess == cudaStat4);
    cudaStat1 = cudaMemcpy(d_A, A, sizeof(double) * lda * n   , cudaMemcpyHostToDevice);
    cudaStat2 = cudaMemcpy(d_B, B, sizeof(double) * ldb * nrhs, cudaMemcpyHostToDevice);
    assert(cudaSuccess == cudaStat1);
    assert(cudaSuccess == cudaStat2);

// step 3: query working space of geqrf and ormqr
    cusolver_status = cusolverDnDgeqrf_bufferSize(
        cudenseH,
        m,
        n,
        d_A,
        lda,
        &lwork);
    assert (cusolver_status == CUSOLVER_STATUS_SUCCESS);

    cudaStat1 = cudaMalloc((void**)&d_work, sizeof(double)*lwork);
    assert(cudaSuccess == cudaStat1);

// step 4: compute QR factorization
    cusolver_status = cusolverDnDgeqrf(
        cudenseH,
        m,
        n,
        d_A,
        lda,
        d_tau,
        d_work,
        lwork,
        devInfo);
    cudaStat1 = cudaDeviceSynchronize();
    assert(CUSOLVER_STATUS_SUCCESS == cusolver_status);
    assert(cudaSuccess == cudaStat1);

    // check if QR is good or not
    cudaStat1 = cudaMemcpy(&info_gpu, devInfo, sizeof(int), cudaMemcpyDeviceToHost);
    assert(cudaSuccess == cudaStat1);

    printf("after geqrf: info_gpu = %d\n", info_gpu);
    assert(0 == info_gpu);

// step 5: compute Q^T*B
    cusolver_status= cusolverDnDormqr(
        cudenseH,
        CUBLAS_SIDE_LEFT,
        CUBLAS_OP_T,
        m,
        nrhs,
        n,
        d_A,
        lda,
        d_tau,
        d_B,
        ldb,
        d_work,
        lwork,
        devInfo);
    cudaStat1 = cudaDeviceSynchronize();
    assert(CUSOLVER_STATUS_SUCCESS == cusolver_status);
    assert(cudaSuccess == cudaStat1);

    // check if QR is good or not
    cudaStat1 = cudaMemcpy(&info_gpu, devInfo, sizeof(int), cudaMemcpyDeviceToHost);
    assert(cudaSuccess == cudaStat1);

    printf("after ormqr: info_gpu = %d\n", info_gpu);
    assert(0 == info_gpu);


// step 6: compute x = R \ Q^T*B

    cublas_status = cublasDtrsm(
         cublasH,
         CUBLAS_SIDE_LEFT,
         CUBLAS_FILL_MODE_UPPER,
         CUBLAS_OP_N,
         CUBLAS_DIAG_NON_UNIT,
         n,
         nrhs,
         &one,
         d_A,
         lda,
         d_B,
         ldb);
    cudaStat1 = cudaDeviceSynchronize();
    assert(CUBLAS_STATUS_SUCCESS == cublas_status);
    assert(cudaSuccess == cudaStat1);

    cudaStat1 = cudaMemcpy(XC, d_B, sizeof(double)*ldb*nrhs, cudaMemcpyDeviceToHost);
    assert(cudaSuccess == cudaStat1);

    printf("X = (matlab base-1)\n");
    printMatrix(n, nrhs, XC, ldb, "X");

// free resources
    if (d_A    ) cudaFree(d_A);
    if (d_tau  ) cudaFree(d_tau);
    if (d_B    ) cudaFree(d_B);
    if (devInfo) cudaFree(devInfo);
    if (d_work ) cudaFree(d_work);


    if (cublasH ) cublasDestroy(cublasH);
    if (cudenseH) cusolverDnDestroy(cudenseH);

    cudaDeviceReset();

    return 0; 
}
