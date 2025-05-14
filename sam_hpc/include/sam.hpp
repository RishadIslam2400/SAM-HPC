#pragma once

#include "CSRMatrix.hpp"
#include "sparsityPattern.hpp"
#include "householderQR.hpp"
#include "mgsQR.hpp"
#include "eigenQRSolve.hpp"

#include "launchThreads.hpp"

// @todo: add post filtration

template <typename T, typename SparsityPatternType>
class SparseApproximateMap {
private:
    const CSRMatrix<T> *targetMatrix;
    const CSRMatrix<T> *sourceMatrix;
    CSRMatrix<T> *mappingMatrix;
    SparsityPattern<T, SparsityPatternType> *sparsityPattern;

    // Helper functions
    
public:
    SparseApproximateMap(const CSRMatrix<T> &targetMatrix, const CSRMatrix<T> &sourceMatrix, const SparsityPattern<T, SparsityPatternType> &pattern);
    SparseApproximateMap(const CSRMatrix<T> &targetMatrix, const CSRMatrix<T> &sourceMatrix, const SparsityPatternType &patternType);

    SparseApproximateMap() = delete;
    SparseApproximateMap(const SparseApproximateMap &other) = delete;
    SparseApproximateMap &operator=(const SparseApproximateMap &other) = delete;
    SparseApproximateMap(SparseApproximateMap &&other) = delete;
    SparseApproximateMap &operator=(SparseApproximateMap &&other) = delete;

    void computeMap();
    const CSRMatrix<T> *getMap() const;
    const CSRMatrix<int> *getPattern() const;

    ~SparseApproximateMap();
};

template <typename T, typename SparsityPatternType>
SparseApproximateMap<T, SparsityPatternType>::SparseApproximateMap(const CSRMatrix<T> &targetMatrix,
                                                                   const CSRMatrix<T> &sourceMatrix,
                                                                   const SparsityPattern<T, SparsityPatternType> &pattern) {
    assert((targetMatrix.row_num == sourceMatrix.row_num && targetMatrix.col_num == sourceMatrix.col_num) && "Target and source matrices must have the same dimensions.");
    assert((targetMatrix.row_num > 0 && targetMatrix.col_num > 0) && "Target and source matrices must have non-zero dimensions.");
    this->targetMatrix = &targetMatrix;
    this->sourceMatrix = &sourceMatrix;
    this->mappingMatrix = nullptr;
    this->sparsityPattern = new SparsityPattern<T, SparsityPatternType>(pattern);
}

template <typename T, typename SparsityPatternType>
SparseApproximateMap<T, SparsityPatternType>::SparseApproximateMap(const CSRMatrix<T> &targetMatrix, const CSRMatrix<T> &sourceMatrix, const SparsityPatternType &patternType) {
    assert((targetMatrix.row_num == sourceMatrix.row_num && targetMatrix.col_num == sourceMatrix.col_num) && "Target and source matrices must have the same dimensions.");
    assert((targetMatrix.row_num > 0 && targetMatrix.col_num > 0) && "Target and source matrices must have non-zero dimensions.");
    this->targetMatrix = &targetMatrix;
    this->sourceMatrix = &sourceMatrix;
    this->mappingMatrix = nullptr;
    this->sparsityPattern = new SparsityPattern<T, SparsityPatternType>(sourceMatrix, patternType);
}

template <typename T, typename SparsityPatternType>
void SparseApproximateMap<T, SparsityPatternType>::computeMap() {
    // Get the computed sparsity pattern matrix
    sparsityPattern->computePattern();
    const CSRMatrix<int> *pattern = sparsityPattern->getPattern();

    // Initialize the mapping matrix
    mappingMatrix = new CSRMatrix<T>();
    mappingMatrix->row_num = targetMatrix->row_num;
    mappingMatrix->col_num = targetMatrix->col_num;
    mappingMatrix->nnz = pattern->nnz;
    mappingMatrix->row_pointers = new std::vector<size_t>(*(pattern->row_pointers));
    mappingMatrix->col_indices = new std::vector<size_t>(*(pattern->col_indices));
    mappingMatrix->vals = new std::vector<T>(mappingMatrix->nnz, 0);
    // Todo: change span to pointers
    // Start SAM computation
    if constexpr (SEQUENTIAL) {
        std::vector<size_t> J;
        std::vector<T> rhs;
        std::vector<std::vector<T>> submatrix;
        for (size_t i = 0; i < mappingMatrix->row_num; ++i) {
            // Compute the submatrix indices for each row
            const size_t rowStart = (*(pattern->row_pointers))[i];
            const size_t rowEnd = (*(pattern->row_pointers))[i + 1];
            size_t iSize = rowEnd - rowStart; // For I
            J.clear();
            std::vector<int> marker(mappingMatrix->col_num, -1);

            for (size_t j = rowStart; j < rowEnd; ++j) {
                const size_t colIdx = (*(pattern->col_indices))[j];

                for (size_t colIdxRowStart = (*(pattern->row_pointers))[colIdx]; colIdxRowStart < (*(pattern->row_pointers))[colIdx + 1]; ++colIdxRowStart) {
                    const size_t submatrixColIdx = (*(pattern->col_indices))[colIdxRowStart];
                    if (marker[submatrixColIdx] < 0) {
                        marker[submatrixColIdx] = 1;
                        J.push_back(submatrixColIdx);
                    }
                }
            }

            std::sort(J.begin(), J.end());

            // Compute the RHS vector from the target matrix
            rhs.assign(J.size(), 0);
            size_t targetRowStart = (*(targetMatrix->row_pointers))[i];
            const size_t targetRowEnd = (*(targetMatrix->row_pointers))[i + 1];
            for (size_t j = 0; j < J.size() && targetRowStart < targetRowEnd;) {
                if (J[j] == (*(targetMatrix->col_indices))[targetRowStart]) {
                    rhs[j] = (*(targetMatrix->vals))[targetRowStart];
                    ++targetRowStart;
                    ++j;
                    continue;
                } else if (J[j] > (*(targetMatrix->col_indices))[targetRowStart]) {
                    ++targetRowStart;
                    continue;
                }
                ++j;
            }

            // Compute the submatrix
            submatrix.assign(iSize, std::vector<T>(J.size(), 0));
            for (size_t submatrixRowIdx = 0, j = rowStart; j < rowEnd; ++submatrixRowIdx, ++j) {
                const size_t colIdx = (*(pattern->col_indices))[j];

                size_t sourceRowStart = (*(sourceMatrix->row_pointers))[colIdx];
                const size_t sourceRowEnd = (*(sourceMatrix->row_pointers))[colIdx + 1];
                for (size_t submatrixColIdx = 0; submatrixColIdx < J.size() && sourceRowStart < sourceRowEnd;) {
                    if (J[submatrixColIdx] == (*(sourceMatrix->col_indices))[sourceRowStart]) {
                        // Column major submatrix - striding by row size
                        submatrix[submatrixRowIdx][submatrixColIdx] = (*(sourceMatrix->vals))[sourceRowStart];
                        ++sourceRowStart;
                        ++submatrixColIdx;
                        continue;
                    } else if (J[submatrixColIdx] > (*(sourceMatrix->col_indices))[sourceRowStart]) {
                        ++sourceRowStart;
                        continue;
                    }
                    ++submatrixColIdx;
                }
            }

            std::span<T> x_subspan(mappingMatrix->vals->data() + rowStart, iSize);
            // Solve the submatrix
            // mgsQRSolve(submatrix, rhs, x.subspan(rowStart, iSize), iSize, J.size());
            householderQRSolve(submatrix, rhs, x_subspan, iSize, J.size());
            // eigenQRSolve(submatrix, rhs, x.subspan(rowStart, iSize), iSize, J.size());
        }
    } else {
        auto samComputation = [&](size_t start, size_t end) {
            std::vector<size_t> J;
            std::vector<T> rhs;
            std::vector<std::vector<T>> submatrix;
            for (size_t i = start; i < end; ++i) {
                // Compute the submatrix indices for each row
                const size_t rowStart = (*(pattern->row_pointers))[i];
                const size_t rowEnd = (*(pattern->row_pointers))[i + 1];
                size_t iSize = rowEnd - rowStart; // For I
                J.clear();
                std::vector<int> marker(mappingMatrix->col_num, -1);

                for (size_t j = rowStart; j < rowEnd; ++j) {
                    const size_t colIdx = (*(pattern->col_indices))[j];

                    for (size_t colIdxRowStart = (*(pattern->row_pointers))[colIdx]; colIdxRowStart < (*(pattern->row_pointers))[colIdx + 1]; ++colIdxRowStart) {
                        const size_t submatrixColIdx = (*(pattern->col_indices))[colIdxRowStart];
                        if (marker[submatrixColIdx] < 0) {
                            marker[submatrixColIdx] = 1;
                            J.push_back(submatrixColIdx);
                        }
                    }
                }

                std::sort(J.begin(), J.end());

                // Compute the RHS vector from the target matrix
                rhs.assign(J.size(), 0);
                size_t targetRowStart = (*(targetMatrix->row_pointers))[i];
                const size_t targetRowEnd = (*(targetMatrix->row_pointers))[i + 1];
                for (size_t j = 0; j < J.size() && targetRowStart < targetRowEnd;) {
                    if (J[j] == (*(targetMatrix->col_indices))[targetRowStart]) {
                        rhs[j] = (*(targetMatrix->vals))[targetRowStart];
                        ++targetRowStart;
                        ++j;
                        continue;
                    } else if (J[j] > (*(targetMatrix->col_indices))[targetRowStart]) {
                        ++targetRowStart;
                        continue;
                    }
                    ++j;
                }

                // Compute the submatrix
                submatrix.assign(iSize, std::vector<T>(J.size(), 0));
                for (size_t submatrixRowIdx = 0, j = rowStart; j < rowEnd; ++submatrixRowIdx, ++j) {
                    const size_t colIdx = (*(pattern->col_indices))[j];

                    size_t sourceRowStart = (*(sourceMatrix->row_pointers))[colIdx];
                    const size_t sourceRowEnd = (*(sourceMatrix->row_pointers))[colIdx + 1];
                    for (size_t submatrixColIdx = 0; submatrixColIdx < J.size() && sourceRowStart < sourceRowEnd;) {
                        if (J[submatrixColIdx] == (*(sourceMatrix->col_indices))[sourceRowStart]) {
                            // Column major submatrix - striding by row size
                            submatrix[submatrixRowIdx][submatrixColIdx] = (*(sourceMatrix->vals))[sourceRowStart];
                            ++sourceRowStart;
                            ++submatrixColIdx;
                            continue;
                        } else if (J[submatrixColIdx] > (*(sourceMatrix->col_indices))[sourceRowStart]) {
                            ++sourceRowStart;
                            continue;
                        }
                        ++submatrixColIdx;
                    }
                }

                std::span<T> x_subspan(mappingMatrix->vals->data() + rowStart, iSize);
                // Solve the submatrix
                // mgsQRSolve(submatrix, rhs, x_subspan, iSize, J.size());
                householderQRSolve(submatrix, rhs, x_subspan, iSize, J.size());
                // eigenQRSolve(submatrix, rhs, x_subspan, iSize, J.size());
            }
        };

        launchThreads(mappingMatrix->row_num, samComputation);
    }
}

template <typename T, typename SparsityPatternType>
const CSRMatrix<T> *SparseApproximateMap<T, SparsityPatternType>::getMap() const {
    assert(mappingMatrix != nullptr && "Compute the mapping matrix first.");
    return mappingMatrix;
}

template <typename T, typename SparsityPatternType>
const CSRMatrix<int> *SparseApproximateMap<T, SparsityPatternType>::getPattern() const {
    assert(sparsityPattern != nullptr && "Compute the sparsity pattern first.");
    return sparsityPattern->getPattern();
}

template <typename T, typename SparsityPatternType>
SparseApproximateMap<T, SparsityPatternType>::~SparseApproximateMap() {
    
    delete sparsityPattern;
    
    if (mappingMatrix != nullptr)
        delete mappingMatrix;
}
