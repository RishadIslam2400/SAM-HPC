#pragma once

#include "CSRMatrix.hpp"
#include "extendPattern.hpp"
#include "launchThreads.hpp"

#include <queue>

struct SimplePattern {};
struct GlobalThresholdPattern { double globalThreshold; };
struct ColumnThresholdPattern { double columnThreshold; };
struct FixedNNZPattern { size_t fixedNNZ; };
struct combinedThresholdPattern { double thresh; };

// @todo: perform set union between the computed sparsity pattern and the target matrix patterns
// so that while computing the map we dont skip any non-zero entries from the target matrix
// @todo: get rid of private
template <typename T, typename PatternType>
class SparsityPattern {
public:
    SparsityPattern() = delete;

    SparsityPattern(const CSRMatrix<T> &originalMatrix, const PatternType &type, int level = 2)
        : m_originalMatrix(originalMatrix), m_type(type), m_level(level), m_pattern(nullptr) {}
    
    SparsityPattern(const SparsityPattern &other) : m_originalMatrix(other.m_originalMatrix), m_type(other.m_type), m_level(other.m_level) {
        if (other.m_pattern)
            m_pattern = std::make_unique<CSRMatrix<int>>(*other.m_pattern);
    }

    SparsityPattern &operator=(const SparsityPattern &other) {
        if (this != &other) {
            m_originalMatrix = other.m_originalMatrix;
            m_type = other.m_type;
            m_level = other.m_level;
            m_pattern.release();
            if (other.m_pattern)
                m_pattern = std::make_unique<CSRMatrix<int>>(*other.m_pattern);
        }

        return *this;
    }

    SparsityPattern(SparsityPattern&&) = default;
    SparsityPattern& operator=(SparsityPattern&&) = default;

    template <typename X, typename Type>
    friend bool operator==(const SparsityPattern<X, Type> &lhs, const SparsityPattern<X, Type> &rhs);

    template <typename X, typename Type>
    friend bool operator!=(const SparsityPattern<X, Type> &lhs, const SparsityPattern<X, Type> &rhs);

    template <typename X, typename Type>
    friend std::ostream &operator<<(std::ostream &os, const SparsityPattern<X, Type> &p);

    void computePattern() {
        if (m_pattern) return; // Pattern already computed

        if constexpr (std::is_same_v<PatternType, SimplePattern>)
            computeSimplePattern();
        else if constexpr (std::is_same_v<PatternType, GlobalThresholdPattern>)
            computeGlobalThresholdPattern(m_type.globalThreshold);
        else if constexpr (std::is_same_v<PatternType, ColumnThresholdPattern>)
            computeColumnThresholdPattern(m_type.columnThreshold);
        else if constexpr (std::is_same_v<PatternType, FixedNNZPattern>)
            computeFixedNNZPattern(m_type.fixedNNZ);
        else if constexpr (std::is_same_v<PatternType, combinedThresholdPattern>)
            computeCombinedPattern(m_type.thresh);
        else
            static_assert(std::is_void_v<T>, "Unsupported pattern type");
    }

    const CSRMatrix<int> *getPattern() const {
        assert(m_pattern != nullptr && "Compute the sparsity pattern first.");
        return m_pattern.get();
    }
    
    size_t getNNZ() const {
        return m_pattern ? m_pattern->m_nnz : 0;
    }

private:
    const CSRMatrix<T> &m_originalMatrix;      // source matrix
    PatternType m_type;                        // sparsification technique
    int m_level;                               // level of pattern extension
    std::unique_ptr<CSRMatrix<int>> m_pattern; // computed pattern

    // ================ Sparsity Pattern Computation ================
    void computeSimplePattern() {
        std::vector<int> patternValues(m_originalMatrix.m_nnz, 1);
        m_pattern = std::make_unique<CSRMatrix<int>>(
            m_originalMatrix.m_rows,
            m_originalMatrix.m_cols,
            patternValues,
            m_originalMatrix.m_row_pointers,
            m_originalMatrix.m_col_indices
        );

        // sam::extend_pattern(*m_pattern, m_level);
    }

    void buildPatternFromFilteredCols(std::vector<std::vector<size_t>>& filtered_cols) {
        m_pattern = std::make_unique<CSRMatrix<int>>();
        m_pattern->m_rows = m_originalMatrix.m_rows;
        m_pattern->m_cols = m_originalMatrix.m_cols;
        m_pattern->m_row_pointers.resize(m_pattern->m_rows + 1, 0);

        // Compute the row pointers
        for (size_t i = 0; i < m_pattern->m_rows; ++i) {
            m_pattern->m_row_pointers[i + 1] = m_pattern->m_row_pointers[i] + filtered_cols[i].size();
        }

        m_pattern->m_nnz = m_pattern->m_row_pointers.back();
        m_pattern->m_col_indices.resize(m_pattern->m_nnz);
        m_pattern->m_vals.resize(m_pattern->m_nnz, 1);

        tbb::parallel_for(tbb::blocked_range<size_t>(0, m_pattern->m_rows),
            [&](const tbb::blocked_range<size_t>& r) {
                for (size_t i = r.begin(); i < r.end(); ++i) {
                    if (!filtered_cols[i].empty()) {
                        size_t *dest = m_pattern->m_col_indices.data() + m_pattern->m_row_pointers[i];
                        std::copy(filtered_cols[i].begin(), filtered_cols[i].end(), dest);
                    }
                 }
            }
        );

        sam::extend_pattern(*m_pattern, m_level);
    }

    void computeGlobalThresholdPattern(const double globalThreshold) {
        const std::vector<T> diag = diagonal<diagonalType::forScaling>(m_originalMatrix);
        std::vector<T> scaledValues = m_originalMatrix.m_vals;
        diagonalScaling(scaledValues, diag);

        std::vector<std::vector<size_t>> filtered_cols(m_originalMatrix.m_rows);

        // Count number of non-zero elements in each row
        tbb::parallel_for(tbb::blocked_range<size_t>(0, m_originalMatrix.m_rows),
            [&](const tbb::blocked_range<size_t>& r) {
                for (size_t i = r.begin(); i < r.end(); ++i) {
                    bool diagonal_found = false;

                    const size_t rowStart = m_originalMatrix.m_row_pointers[i];
                    const size_t rowEnd = m_originalMatrix.m_row_pointers[i + 1];
                    filtered_cols[i].reserve(rowEnd - rowStart);
                    for (size_t j = rowStart; j < rowEnd; ++j) {
                        const size_t colIdx = m_originalMatrix.m_col_indices[j];
                        if (colIdx == i) {
                            diagonal_found = true;
                            filtered_cols[i].push_back(colIdx);
                        } else if (std::abs(scaledValues[j]) > globalThreshold) {
                            filtered_cols[i].push_back(colIdx);
                        }
                    }

                    if (!diagonal_found) {
                        filtered_cols[i].push_back(i);
                    }

                    std::sort(filtered_cols[i].begin(), filtered_cols[i].end());
                }
            }
        );

        buildPatternFromFilteredCols(filtered_cols);
    }

    void computeColumnThresholdPattern(const double tau) {
        std::vector<std::vector<size_t>> filtered_cols(m_originalMatrix.m_rows); 

        // Count number of non-zero elements in each row
        tbb::parallel_for(tbb::blocked_range<size_t>(0, m_originalMatrix.m_rows),
            [&](const tbb::blocked_range<size_t>& r) {
                for (size_t i = r.begin(); i < r.end(); ++i) {
                    const size_t rowStart = m_originalMatrix.m_row_pointers[i];
                    const size_t rowEnd = m_originalMatrix.m_row_pointers[i + 1];

                    // Find max absolute value in the row
                    T maxVal = 0;
                    for (size_t j = rowStart; j < rowEnd; ++j) {
                        maxVal = std::max(maxVal, std::abs(m_originalMatrix.m_vals[j]));
                    }
                    const T threshold = (1 - tau) * maxVal;

                    // Count non-zeros using the just-computed threshold
                    bool diagonal_found = false;
                    filtered_cols[i].reserve(rowEnd - rowStart);
                    for (size_t j = rowStart; j < rowEnd; ++j) {
                        const size_t colIdx = m_originalMatrix.m_col_indices[j];
                        if (colIdx == i) {
                            diagonal_found = true;
                            filtered_cols[i].push_back(colIdx);
                        } else if (std::abs(m_originalMatrix.m_vals[j]) > threshold) {
                            filtered_cols[i].push_back(colIdx);
                        }
                    }

                    if (!diagonal_found) {
                        filtered_cols[i].push_back(i);
                    }

                    std::sort(filtered_cols[i].begin(), filtered_cols[i].end());
                }
            }
        );

        buildPatternFromFilteredCols(filtered_cols);
    }

    void computeFixedNNZPattern(const size_t lfil) {
        std::vector<std::vector<size_t>> filtered_cols(m_originalMatrix.m_rows);

        // Count number of non-zero elements in each row
        tbb::parallel_for(tbb::blocked_range<size_t>(0, m_originalMatrix.m_rows), [&](const tbb::blocked_range<size_t>& r) {
            for (size_t i = r.begin(); i < r.end(); ++i) {
                const size_t rowStart = m_originalMatrix.m_row_pointers[i];
                const size_t rowEnd = m_originalMatrix.m_row_pointers[i + 1];
                const size_t nnz = rowEnd - rowStart;

                if (nnz == 0) continue;

                std::vector<std::pair<T, size_t>> col_entries;
                col_entries.reserve(nnz);

                for (size_t j = rowStart; j < rowEnd; ++j) {
                    col_entries.emplace_back(std::abs(m_originalMatrix.m_vals[j]), m_originalMatrix.m_col_indices[j]);
                }

                size_t count = std::min(lfil, nnz);
                if (count == 0 && nnz > 0) continue;

                // Partitions the column indices, count is the pivot
                std::nth_element(col_entries.begin(), col_entries.begin() + count - 1, col_entries.end(), std::greater<>{});

                filtered_cols[i].reserve(count);
                for (size_t k = 0; k < count; ++k) {
                    filtered_cols[i].push_back(col_entries[k].second);
                }

                std::sort(filtered_cols[i].begin(), filtered_cols[i].end());
            }
        });

        buildPatternFromFilteredCols(filtered_cols);
    }
    void computeCombinedPattern(const double thresh);

    // =============== Helper Functions ================
    void diagonalScaling(std::vector<T> &values, const std::vector<T> &diagonal) {
        const std::vector<size_t> &rowPointers = m_originalMatrix.m_row_pointers;
        const std::vector<size_t> &colIndices = m_originalMatrix.m_col_indices;

        // Pre multiplying and post multiplying - (D^-1/2 * A * D^-1/2)
        // multiply each diagonal element with corresponding row in the matrix
        // [a_i] = d_i * [a_i]
        // multiply each diagonal element with corresponding column in the matrix
        // Multiplying each row elements with their corresponding diagonal element (same idx)
        // [a_i]^T = d_i * [a_i]^T
        // @todo: this loop can be done in parallel
        tbb::parallel_for(tbb::blocked_range<size_t>(0, m_originalMatrix.m_rows), [&](const tbb::blocked_range<size_t> &r) {
            for (size_t i = r.begin(); i < r.end(); ++i) {
                const size_t rowStart = rowPointers[i];
                const size_t rowEnd = rowPointers[i + 1];
                for (size_t j = rowStart; j < rowEnd; ++j) {
                    size_t idx = colIndices[j];
                    values[j] *= diagonal[i] * diagonal[idx]; // diagonal[idx] is a random memory access
                }
            }
        });
    }
};

template <typename T, typename PatternType>
void SparsityPattern<T, PatternType>::computeCombinedPattern(const double thresh) {
    // Perform diagonal scaling
    // compute global threshold from local thresholds
    // filter rows with threshold values
}

template <typename X, typename PatternType>
bool operator==(const SparsityPattern<X, PatternType> &lhs, const SparsityPattern<X, PatternType> &rhs) {
    return ((*(lhs.originalMatrix) == *(rhs.originalMatrix)) &&
            ((lhs.pattern == nullptr && rhs.pattern == nullptr) ||
             (lhs.pattern != nullptr && rhs.pattern != nullptr && *(lhs.pattern) == *(rhs.pattern))));
}

template <typename X, typename PatternType>
bool operator!=(const SparsityPattern<X, PatternType> &lhs, const SparsityPattern<X, PatternType> &rhs) {
    return !(lhs == rhs);
}

template <typename X, typename Type>
std::ostream &operator<<(std::ostream &os, const SparsityPattern<X, Type> &p) {
    os << "Sparsity Pattern: " << std::endl;
    os << *(p.pattern);
    return os;
}