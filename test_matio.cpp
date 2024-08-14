#include <vector>
#include <iostream>
#include <mat.h>

int main(int argc, char **argv) 
{
    MATFile *matfile = matOpen("matrix_2.mat", "r");
    if (!matfile) {
        std::cerr << "Error opening file" << std::endl;
        return -1;
    }

    mxArray *jac = matGetVariable(matfile, "Jac");
    if (!jac) {
        std::cerr << "Error getting variable" << std::endl;
        matClose(matfile);
        return -1;
    }

    if (!mxIsSparse(jac))
    {
        std::cerr << "The variable is not a sparse matrix!" << std::endl;
        mxDestroyArray(jac);
        matClose(matfile);
        return -1;
    }

    // Get pointers to the sparse matrix data
    double *pr = mxGetPr(jac);    // Non-zero elements
    mwIndex *ir = mxGetIr(jac);   // Row indices
    mwIndex *jc = mxGetJc(jac);   // Column pointers
    mwSize numRows = mxGetM(jac); // Number of rows
    mwSize numCols = mxGetN(jac); // Number of columns
    mwSize nnz = jc[numCols];     // Number of non-zero elements

    std::cout << "Number of rows: " << numRows << std::endl;
    std::cout << "Number of columns: " << numCols << std::endl;
    std::cout << "Number of non-zero elements: " << nnz << std::endl;

    // Cleanup
    mxDestroyArray(jac);
    matClose(matfile);
}