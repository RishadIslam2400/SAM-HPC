#include "CSRMatrix.hpp"
#include "sparsityPattern.hpp"

template <typename T>
class sparseApproximateMap
{
private:
    const SparseMatrix::CSRMatrix<T> *targetMatrix;
    const SparseMatrix::CSRMatrix<T> *sourceMatrix;
    SparseMatrix::CSRMatrix<T> *mappingMatrix;
    SparsityPattern<T> *sparsityPattern;    

public:
    sparseApproximateMap() = delete;
    sparseApproximateMap(const SparseMatrix::CSRMatrix<T> &targetMatrix, const SparseMatrix::CSRMatrix<T> &sourceMatrix);

    ~sparseApproximateMap();
};

template <typename T>
sparseApproximateMap<T>::sparseApproximateMap(const SparseMatrix::CSRMatrix<T> &targetMatrix, const SparseMatrix::CSRMatrix<T> &sourceMatrix)
    : targetMatrix(&targetMatrix), sourceMatrix(&sourceMatrix)
{
    mappingMatrix = nullptr;
    sparsityPattern = nullptr;
}

template <typename T>
sparseApproximateMap<T>::~sparseApproximateMap()
{
    delete mappingMatrix;
    delete sparsityPattern;
}
