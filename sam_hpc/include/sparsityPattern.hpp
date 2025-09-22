#pragma once

#include "CSRMatrix.hpp"
#include "extendPattern.hpp"
#include "launchThreads.hpp"

#include <queue>

struct SimplePattern {};
struct GlobalThresholdPattern { double globalThreshold; };
struct ColumnThresholdPattern { double columnThreshold; };
struct FixedNNZPattern { size_t fixedNNZ; };
struct combinedThresholdPattern { 
    double global_thresh;
    double column_thresh;
};

// @todo: perform set union between the computed sparsity pattern and the target matrix patterns
// so that while computing the map we dont skip any non-zero entries from the target matrix
// @todo: get rid of private
template <typename T, typename PatternType>
class SparsityPattern {
public:
    SparsityPattern() = delete;

    SparsityPattern(const CSRMatrix<T> &originalMatrix, const CSRMatrix<T> &targetMatrix, const PatternType &type, int level = 2)
        : m_originalMatrix(originalMatrix), m_targetMatrix(targetMatrix), m_type(type), m_level(level), m_pattern(nullptr) {}
    
    SparsityPattern(const SparsityPattern &other) = delete;
    SparsityPattern &operator=(const SparsityPattern &other) = delete;

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
            computeCombinedPattern(m_type.global_thresh, m_type.column_thresh);
        else
            static_assert(!std::is_same_v<T,T>, "Unsupported pattern type");
    }

    const CSRMatrix<int> *getPattern() const {
        assert(m_pattern != nullptr && "Compute the sparsity pattern first.");
        return m_pattern.get();
    }
    
    size_t getNNZ() const {
        return m_pattern ? m_pattern->m_nnz : 0;
    }

    // For pattern union after extension
    std::shared_ptr<CSRMatrix<int>> pattern_union() {
        // Integrate the target matrix into the pattern and produce the final pattern
        assert(m_pattern && "Initial pattern must be computed before union.");
        
        std::shared_ptr<CSRMatrix<int>> final_pattern = std::make_shared<CSRMatrix<int>>();
        final_pattern->m_rows = m_pattern->m_rows;
        final_pattern->m_cols = m_pattern->m_cols;
        final_pattern->m_row_pointers.resize(final_pattern->m_rows + 1, 0);

        tbb::enumerable_thread_specific<std::vector<int>> local_markers(final_pattern->m_cols, -1);
        tbb::parallel_for(tbb::blocked_range<size_t>(0, final_pattern->m_rows), [&](const tbb::blocked_range<size_t>& r) {
            std::vector<int>& marker = local_markers.local();
            for (size_t i = r.begin(); i < r.end(); ++i) {
                size_t cols = 0;

                for (size_t j = m_pattern->m_row_pointers[i], e = m_pattern->m_row_pointers[i + 1]; j < e; ++j) {
                    const size_t colIdx = m_pattern->m_col_indices[j];
                    if (marker[colIdx] != static_cast<int>(i)) {
                        marker[colIdx] = static_cast<int>(i);
                        cols++;
                    }
                }

                for (size_t j = m_targetMatrix.m_row_pointers[i], e = m_targetMatrix.m_row_pointers[i + 1]; j < e; ++j) {
                    const size_t colIdx = m_targetMatrix.m_col_indices[j];
                    if (marker[colIdx] != static_cast<int>(i)) {
                        marker[colIdx] = static_cast<int>(i);
                        cols++;
                    }
                }

                final_pattern->m_row_pointers[i + 1] = cols;
            }
        });

        final_pattern->m_nnz = final_pattern->scanRowSize();
        final_pattern->m_col_indices.resize(final_pattern->m_nnz, 0);
        final_pattern->m_vals.resize(final_pattern->m_nnz, 1);

        local_markers.clear();
        tbb::parallel_for(tbb::blocked_range<size_t>(0, final_pattern->m_rows), [&](const tbb::blocked_range<size_t>& r) {
            std::vector<int>& marker = local_markers.local();
            for (size_t i = r.begin(); i < r.end(); ++i) {
                const size_t rowBeg = final_pattern->m_row_pointers[i];
                size_t rowEnd = rowBeg;

                for (size_t j = m_pattern->m_row_pointers[i], e = m_pattern->m_row_pointers[i + 1]; j < e; ++j) {
                    size_t colIdx = m_pattern->m_col_indices[j];

                    if (marker[colIdx] < static_cast<int>(rowBeg)) {
                        marker[colIdx] = static_cast<int>(rowEnd);
                        final_pattern->m_col_indices[rowEnd] = colIdx;
                        rowEnd++;
                    }
                }

                for (size_t j = m_targetMatrix.m_row_pointers[i], e = m_targetMatrix.m_row_pointers[i + 1]; j < e; ++j) {
                    size_t colIdx = m_targetMatrix.m_col_indices[j];

                    if (marker[colIdx] < static_cast<int>(rowBeg)) {
                        marker[colIdx] = static_cast<int>(rowEnd);
                        final_pattern->m_col_indices[rowEnd] = colIdx;
                        rowEnd++;
                    }
                }

                // Clean up the marker entries used for the current row.
                /* for (size_t k = rowBeg; k < rowEnd; ++k) {
                    marker[final_pattern->m_col_indices[k]] = -1;
                } */

                
                sortRow(final_pattern->m_col_indices.data() + rowBeg, final_pattern->m_vals.data() + rowBeg, static_cast<int>(rowEnd - rowBeg));
            }
        });

        return final_pattern;
    }

private:
    const CSRMatrix<T> &m_originalMatrix;      // source matrix
    const CSRMatrix<T> &m_targetMatrix;        // target matrix
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

    template<typename Func>
    void buildPattern(Func filter) {
        m_pattern = std::make_unique<CSRMatrix<int>>();
        m_pattern->m_rows = m_originalMatrix.m_rows;
        m_pattern->m_cols = m_originalMatrix.m_cols;
        m_pattern->m_row_pointers.resize(m_pattern->m_rows + 1, 0);

        // Count nnz per row
        tbb::parallel_for(tbb::blocked_range<size_t>(0, m_pattern->m_rows), [&](const tbb::blocked_range<size_t> &r) {
            for (size_t i = r.begin(); i < r.end(); ++i) {
                m_pattern->m_row_pointers[i + 1] = filter(i, nullptr); // no write operation
            }
        });

        std::partial_sum(m_pattern->m_row_pointers.begin(), m_pattern->m_row_pointers.end(), m_pattern->m_row_pointers.begin());
        m_pattern->m_nnz = m_pattern->m_row_pointers.back();
        m_pattern->m_col_indices.resize(m_pattern->m_nnz, 0);
        m_pattern->m_vals.resize(m_pattern->m_nnz, 1);

        // Fill the column indices
        tbb::parallel_for(tbb::blocked_range<size_t>(0, m_originalMatrix.m_rows), [&](const tbb::blocked_range<size_t> &r) {
            for (size_t i = r.begin(); i < r.end(); ++i) {
                size_t *dest = m_pattern->m_col_indices.data() + m_pattern->m_row_pointers[i];
                filter(i, dest);
                std::sort(dest, m_pattern->m_col_indices.data() + m_pattern->m_row_pointers[i + 1]);
            }
        });

        sam::extend_pattern(*m_pattern, m_level);
    }

    void computeGlobalThresholdPattern(const double globalThreshold) {
        const std::vector<T> diag = diagonal<diagonalType::forScaling>(m_originalMatrix);
        std::vector<T> scaledValues = m_originalMatrix.m_vals;
        diagonalScaling(scaledValues, diag);

        auto filter = [&](size_t i, size_t* dest) {
            size_t count = 0;
            bool diagonal_found = false;
            for (size_t j = m_originalMatrix.m_row_pointers[i]; j < m_originalMatrix.m_row_pointers[i + 1]; ++j) {
                const size_t colIdx = m_originalMatrix.m_col_indices[j];
                bool keep = false;
                if (colIdx == i) {
                    keep = true;
                    diagonal_found = true;
                } else if (std::abs(scaledValues[j]) > globalThreshold) {
                    keep = true;
                }
                if (keep) {
                    if (dest) dest[count] = colIdx;
                    count++;
                }
            }
            if (!diagonal_found) {
                if (dest) dest[count] = i; // Ensure diagonal is present
                count++;
            }
            return count;
        };

        buildPattern(filter);
    }

    void computeColumnThresholdPattern(const double tau) {
        const std::vector<T> diag = diagonal<diagonalType::forScaling>(m_originalMatrix);
        std::vector<T> scaledValues = m_originalMatrix.m_vals;
        diagonalScaling(scaledValues, diag);

        auto filter = [&](size_t i, size_t* dest) {
            const size_t rowStart = m_originalMatrix.m_row_pointers[i];
            const size_t rowEnd = m_originalMatrix.m_row_pointers[i + 1];

            // Find max absolute value in the row
            T maxVal = 0;
            for (size_t j = rowStart; j < rowEnd; ++j) {
                maxVal = std::max(maxVal, std::abs(scaledValues[j]));
            }
            const T threshold = (1 - tau) * maxVal;


            size_t count = 0;
            bool diagonal_found = false;
            for (size_t j = rowStart; j < rowEnd; ++j) {
                const size_t colIdx = m_originalMatrix.m_col_indices[j];
                bool keep = false;
                if (colIdx == i) {
                    keep = true;
                    diagonal_found = true;
                } else if (std::abs(scaledValues[j]) > threshold) {
                    keep = true;
                }
                if (keep) {
                    if (dest) dest[count] = colIdx;
                    count++;
                }
            }
            if (!diagonal_found) {
                if (dest) dest[count] = i; // Ensure diagonal is present
                count++;
            }
            return count;
        };
        
        buildPattern(filter);
    }

    void computeFixedNNZPattern(const size_t lfil) {
        const std::vector<T> diag = diagonal<diagonalType::forScaling>(m_originalMatrix);
        std::vector<T> scaledValues = m_originalMatrix.m_vals;
        diagonalScaling(scaledValues, diag);

        tbb::enumerable_thread_specific<std::vector<std::pair<T, size_t>>> local_entries;
        auto filter = [&](size_t i, size_t *dest) {
            const size_t rowStart = m_originalMatrix.m_row_pointers[i];
            const size_t rowEnd = m_originalMatrix.m_row_pointers[i + 1];
            const size_t nnz = rowEnd - rowStart;

            size_t count = std::min(lfil, nnz);

            if (dest) {
                auto &col_entries = local_entries.local();
                col_entries.clear();
                col_entries.reserve(nnz);

                for (size_t j = rowStart; j < rowEnd; ++j) {
                    col_entries.emplace_back(std::abs(scaledValues[j]), m_originalMatrix.m_col_indices[j]);
                }

                std::nth_element(col_entries.begin(), col_entries.begin() + count - 1, col_entries.end(), std::greater<>{});

                for (size_t k = 0; k < count; ++k) {
                    dest[k] = col_entries[k].second;
                }
            }

            return count;
        };

        buildPattern(filter);
    }

    void computeCombinedPattern(const double global_thresh, const double column_thresh) {
        const std::vector<T> diag = diagonal<diagonalType::forScaling>(m_originalMatrix);
        std::vector<T> scaledValues = m_originalMatrix.m_vals;
        diagonalScaling(scaledValues, diag);

        auto filter = [&](size_t i, size_t* dest) {
            const size_t rowStart = m_originalMatrix.m_row_pointers[i];
            const size_t rowEnd = m_originalMatrix.m_row_pointers[i + 1];

            // Find max absolute value in the row
            T maxVal = 0;
            for (size_t j = rowStart; j < rowEnd; ++j) {
                maxVal = std::max(maxVal, std::abs(scaledValues[j]));
            }
            const T threshold = (1 - column_thresh) * maxVal;


            size_t count = 0;
            bool diagonal_found = false;
            for (size_t j = rowStart; j < rowEnd; ++j) {
                const size_t colIdx = m_originalMatrix.m_col_indices[j];
                bool keep = false;
                if (colIdx == i) {
                    keep = true;
                    diagonal_found = true;
                } else if (std::abs(scaledValues[j]) > threshold && std::abs(scaledValues[j]) > global_thresh) {
                    keep = true;
                }
                if (keep) {
                    if (dest) dest[count] = colIdx;
                    count++;
                }
            }
            if (!diagonal_found) {
                if (dest) dest[count] = i; // Ensure diagonal is present
                count++;
            }
            return count;
        };

        buildPattern(filter);
    }

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