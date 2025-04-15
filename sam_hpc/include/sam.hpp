#pragma once

#include "CSRMatrix.hpp"
#include "sparsityPattern.hpp"

template <typename T, typename SparsityPatternType>
class SparseApproximateMap
{
private:
    const SparseMatrix::CSRMatrix<T> *targetMatrix;
    const SparseMatrix::CSRMatrix<T> *sourceMatrix;
    SparseMatrix::CSRMatrix<T> *mappingMatrix;
    SparsityPattern<T, SparsityPatternType> *sparsityPattern;

    // Helper functions
    
public:
    SparseApproximateMap(const SparseMatrix::CSRMatrix<T> &targetMatrix, const SparseMatrix::CSRMatrix<T> &sourceMatrix, const SparsityPattern<T, SparsityPatternType> &pattern);
    SparseApproximateMap(const SparseMatrix::CSRMatrix<T> &targetMatrix, const SparseMatrix::CSRMatrix<T> &sourceMatrix, const SparsityPatternType &patternType);

    SparseApproximateMap() = delete;
    SparseApproximateMap(const SparseApproximateMap &other) = delete;
    SparseApproximateMap &operator=(const SparseApproximateMap &other) = delete;
    SparseApproximateMap(SparseApproximateMap &&other) = delete;
    SparseApproximateMap &operator=(SparseApproximateMap &&other) = delete;

    void computeMap();
    const SparseMatrix::CSRMatrix<T> *getMap() const;
    const SparseMatrix::CSRMatrix<int> *getPattern() const;

    ~SparseApproximateMap();
};

template <typename T, typename SparsityPatternType>
SparseApproximateMap<T, SparsityPatternType>::SparseApproximateMap(const SparseMatrix::CSRMatrix<T> &targetMatrix,
                                                                   const SparseMatrix::CSRMatrix<T> &sourceMatrix,
                                                                   const SparsityPattern<T, SparsityPatternType> &pattern)
{
    assert((targetMatrix.row_num == sourceMatrix.row_num && targetMatrix.col_num == sourceMatrix.col_num) && "Target and source matrices must have the same dimensions.");
    assert((targetMatrix.row_num > 0 && targetMatrix.col_num > 0) && "Target and source matrices must have non-zero dimensions.");
    this->targetMatrix = &targetMatrix;
    this->sourceMatrix = &sourceMatrix;
    this->mappingMatrix = nullptr;
    this->sparsityPattern = new SparsityPattern<T, SparsityPatternType>(pattern);
}

template <typename T, typename SparsityPatternType>
SparseApproximateMap<T, SparsityPatternType>::SparseApproximateMap(const SparseMatrix::CSRMatrix<T> &targetMatrix, const SparseMatrix::CSRMatrix<T> &sourceMatrix, const SparsityPatternType &patternType)
{
    assert((targetMatrix.row_num == sourceMatrix.row_num && targetMatrix.col_num == sourceMatrix.col_num) && "Target and source matrices must have the same dimensions.");
    assert((targetMatrix.row_num > 0 && targetMatrix.col_num > 0) && "Target and source matrices must have non-zero dimensions.");
    this->targetMatrix = &targetMatrix;
    this->sourceMatrix = &sourceMatrix;
    this->mappingMatrix = nullptr;
    this->sparsityPattern = new SparsityPattern<T, SparsityPatternType>(sourceMatrix, patternType);
}

template <typename T, typename SparsityPatternType>
void SparseApproximateMap<T, SparsityPatternType>::computeMap()
{
    // Initialize the mapping matrix
    mappingMatrix = new SparseMatrix::CSRMatrix<T>();
    mappingMatrix->row_num = targetMatrix->row_num;
    mappingMatrix->col_num = targetMatrix->col_num;
    mappingMatrix->row_pointers = new std::vector<size_t>(targetMatrix->row_num + 1);

    // Get the computed sparsity pattern matrix
    sparsityPattern->computePattern();
    const SparseMatrix::CSRMatrix<int> *pattern = sparsityPattern->getPattern();

    std::vector<int> marker(mappingMatrix->col_num, -1);
    std::vector<size_t> I, J;
    std::vector<T> rhs;
    std::vector<std::vector<T>> submatrix;
    for (size_t i = 0; i < mappingMatrix->row_num; ++i)
    {
        // Compute the submatrix indices for each row
        const size_t rowStart = (*(pattern->row_pointers))[i];
        const size_t rowEnd = (*(pattern->row_pointers))[i + 1];
        I.assign(pattern->col_indices->begin() + rowStart, pattern->col_indices->begin() + rowEnd);
        J.clear();

        for (size_t j = rowStart; j < rowEnd; ++j)
        {
            const size_t colIdx = (*(pattern->col_indices))[j];

            for (size_t colIdxRowStart = (*(sourceMatrix->row_pointers))[colIdx]; colIdxRowStart < (*(sourceMatrix->row_pointers))[colIdx + 1]; ++colIdxRowStart)
            {
                const size_t submatrixColIdx = (*(sourceMatrix->col_indices))[colIdxRowStart];
                if (marker[submatrixColIdx] < 0)
                {
                    marker[submatrixColIdx] = 1;
                    J.push_back(submatrixColIdx);
                }
            }
        }

        std::sort(J.begin(), J.end());

        // Compute the RHS vector from the target matrix
        rhs.assign(J.size(), 0);
        size_t targetRowStart = (*(targetMatrix->row_pointers))[i];
        size_t targetRowEnd = (*(targetMatrix->row_pointers))[i + 1];
        for (size_t j = 0; j < J.size() && targetRowStart < targetRowEnd;)
        {
            if (J[j] == (*(targetaMatrix->col_indices))[targetRowStart])
            {
                rhs[j] = (*(targetMatrix->vals))[targetRowStart];
                marker[J[j]] = static_cast<int>(j);
                ++targetRowStart;
                ++j;
                continue;
            }
            else if (J[j] < (*(targetMatrix->col_indices))[targetRowStart])
            {
                ++targetRowStart;
                continue;
            }
            ++j;
        }

        // Compute the submatrix
        submatrix.assign(I.size(), std::vector<T>(J.size(), 0));
        for (size_t j = rowStart; j < rowEnd; ++j)
        {
            
        }
    }
}

template <typename T, typename SparsityPatternType>
const SparseMatrix::CSRMatrix<T> *SparseApproximateMap<T, SparsityPatternType>::getMap() const
{
    assert(mappingMatrix != nullptr && "Compute the mapping matrix first.");
    return mappingMatrix;
}

template <typename T, typename SparsityPatternType>
const SparseMatrix::CSRMatrix<int> *SparseApproximateMap<T, SparsityPatternType>::getPattern() const
{
    assert(sparsityPattern != nullptr && "Compute the sparsity pattern first.");
    return sparsityPattern->getPattern();
}

template <typename T, typename SparsityPatternType>
SparseApproximateMap<T, SparsityPatternType>::~SparseApproximateMap()
{
    
    delete sparsityPattern;
    
    if (mappingMatrix != nullptr)
        delete mappingMatrix;
}
