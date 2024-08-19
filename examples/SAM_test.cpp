#include <mat.h>
#include <matrix.h>
#include <iostream>
#include <vector>
#include <tuple>

class csc_matrix
{
public:
    csc_matrix() = default;
    csc_matrix(const std::vector<double>& values, const std::vector<size_t>& rowIndices, const std::vector<size_t>& colPointers, size_t numCols, size_t numRows, size_t nnz)
        : mValues(values), mRowIndices(rowIndices), mColPointers(colPointers), mNumCols(numCols), mNumRows(numRows), mNNZ(nnz) { }

    // Getters
    const std::vector<double>& getValues() const { return mValues; }
    const std::vector<size_t>& getRowIndices() const { return mRowIndices; }
    const std::vector<size_t>& getColPointers() const { return mColPointers; }
    size_t getNumCols() const { return mNumCols; }
    size_t getNumRows() const { return mNumRows; }
    size_t getNNZ() const { return mNNZ; }

    // Setters
    void setValues(const std::vector<double>& values) { mValues = values; }
    void setRowIndices(const std::vector<size_t>& rowIndices) { mRowIndices = rowIndices; }
    void setColPointers(const std::vector<size_t>& colPointers) { mColPointers = colPointers; }
    void setNumCols(size_t numCols) { mNumCols = numCols; }
    void setNumRows(size_t numRows) { mNumRows = numRows; }
    void setNNZ(size_t nnz) { mNNZ = nnz; }

    void printMatrix() const {
        for (size_t col = 0; col < mNumCols; ++col) {
            for (size_t idx = mColPointers[col]; idx < mColPointers[col + 1]; ++idx) {
                std::cout << "(" << mRowIndices[idx] + 1 << ", " << col + 1 << ") = " << mValues[idx] << std::endl;
            }
        }
    }
private:
    std::vector<double> mValues;
    std::vector<size_t> mRowIndices;
    std::vector<size_t> mColPointers;
    size_t mNumCols;
    size_t mNumRows;
    size_t mNNZ;
};

bool read_mat(const char *filename, csc_matrix& sparse_matrix)
{
    // Read the .mat file
    MATFile *matfile = matOpen(filename, "r");
    if (!matfile)
    {
        std::cerr << "Error opening file" << std::endl;
        return false;
    }

    // Read the sparse matrix variable
    mxArray  *a = matGetVariable(matfile, "Jac");
    if (!a) {
        std::cerr << "Error getting variable" << std::endl;
        matClose(matfile);
        return false;
    }

    // Ensure the variable is a sparse matrix
    if (!mxIsSparse(a)) {
        std::cerr << "The variable is not a sparse matrix!" << std::endl;
        mxDestroyArray(a);
        matClose(matfile);
        return false;
    }

    // Get the data
    double *val = mxGetPr(a);                               // Non-zero elements
    mwIndex *rowIndices = mxGetIr(a);                       // Row indices
    mwIndex *colPointers = mxGetJc(a);                      // Column pointers
    size_t numRows = static_cast<size_t>(mxGetM(a));        // Number of rows
    size_t numCols = static_cast<size_t>(mxGetN(a));        // Number of columns
    size_t nnz = static_cast<size_t>(colPointers[numCols]); // Number of non-zero elements

    // Populate the sparse matrix
    std::vector<double> values(val, val + nnz);
    std::vector<size_t> rows(rowIndices, rowIndices + nnz);
    std::vector<size_t> cols(colPointers, colPointers + numCols + 1);
    
    sparse_matrix.setValues(values);
    sparse_matrix.setRowIndices(rows);
    sparse_matrix.setColPointers(cols);
    sparse_matrix.setNumCols(numCols);
    sparse_matrix.setNumRows(numRows);
    sparse_matrix.setNNZ(nnz);

    mxDestroyArray(a);
    matClose(matfile);
    return true;
}

int main()
{
    std::vector<csc_matrix> sequence;
    unsigned int numMatrices = 70;

    csc_matrix A0;
    bool readA0 = read_mat("/home/rishad/SAM-HPC/data/matrix_2.mat", A0);

    if (!readA0) {
        std::cerr << "Error reading matrix file" << std::endl;
        return -1;
    }

    A0.printMatrix();

    return 0;    
}