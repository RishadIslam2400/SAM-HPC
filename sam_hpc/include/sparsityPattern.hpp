#pragma once

#include "CSRMatrix.hpp"
#include "extendPattern.hpp"
#include "helpers.hpp"
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
private:
    CSRMatrix<int> *pattern;
    const CSRMatrix<T> *originalMatrix;
    PatternType type;
    int level;

    // ================ Sparsity Pattern Computation ================
    void computeSimplePattern();
    void computeGlobalThresholdPattern(const double globalThreshold);
    void computeColumnThresholdPattern(const double tau);
    void computeFixedNNZPattern(const size_t lfil);
    void computeCombinedPattern(const double thresh);

    // =============== Helper Functions ================
    void diagonalScaling(std::vector<T> &values, const std::vector<T> &diagonal);

public:
    SparsityPattern() = delete;
    SparsityPattern(const CSRMatrix<T> &originalMatrix, const PatternType &type);
    SparsityPattern(const CSRMatrix<T> &originalMatrix, const PatternType &type, const int level);
    SparsityPattern(const SparsityPattern &other);
    SparsityPattern &operator=(const SparsityPattern &other);
    SparsityPattern(const SparsityPattern &&other);
    SparsityPattern &operator=(const SparsityPattern &&other);

    template <typename X, typename Type>
    friend bool operator==(const SparsityPattern<X, Type> &lhs, const SparsityPattern<X, Type> &rhs);

    template <typename X, typename Type>
    friend bool operator!=(const SparsityPattern<X, Type> &lhs, const SparsityPattern<X, Type> &rhs);

    template <typename X, typename Type>
    friend std::ostream &operator<<(std::ostream &os, const SparsityPattern<X, Type> &p);

    void computePattern();
    const CSRMatrix<int> *getPattern() const;
    size_t getNNZ() const;

    ~SparsityPattern();
};

// ================ Constructor/Assignment ================
template <typename T, typename PatternType>
SparsityPattern<T, PatternType>::SparsityPattern(const CSRMatrix<T> &originalMatrix, const PatternType &type) {
    this->originalMatrix = &originalMatrix;
    this->type = type;
    this->pattern = nullptr;
    this->level = 2;
}

template <typename T, typename PatternType>
SparsityPattern<T, PatternType>::SparsityPattern(const CSRMatrix<T> &originalMatrix, const PatternType &type, const int level) {
    this->originalMatrix = &originalMatrix;
    this->type = type;
    this->pattern = nullptr;
    this->level = level;
}

template <typename T, typename PatternType>
SparsityPattern<T, PatternType>::SparsityPattern(const SparsityPattern &other) {
    originalMatrix = other.originalMatrix;
    type = other.type;
    level = other.level;
    if (other.pattern != nullptr)
        pattern = new CSRMatrix<int>(*other.pattern);
    else
        pattern = nullptr;
}

template <typename T, typename PatternType>
SparsityPattern<T, PatternType> &SparsityPattern<T, PatternType>::operator=(const SparsityPattern &other) {
    if (this != &other) {
        originalMatrix = other.originalMatrix;
        type = other.type;
        level = other.level;
        delete pattern; // Free existing pattern
        if (other.pattern != nullptr)
            pattern = new CSRMatrix<int>(*other.pattern);
        else
            pattern = nullptr;
    }
    return *this;
}

template <typename T, typename PatternType>
SparsityPattern<T, PatternType>::SparsityPattern(const SparsityPattern &&other) {
    originalMatrix = other.originalMatrix;
    type = other.type;
    level = other.level;
    pattern = other.pattern;

    // Set the moved-from object to nullptr
    other.pattern = nullptr;
}

template <typename T, typename PatternType>
SparsityPattern<T, PatternType> &SparsityPattern<T, PatternType>::operator=(const SparsityPattern &&other) {
    if (this != &other) {
        originalMatrix = other.originalMatrix;
        type = other.type;
        level = other.level;

        delete pattern; // Free existing pattern
        pattern = other.pattern;

        // Set the moved-from object to nullptr
        other.pattern = nullptr;
    }
    return *this;
}

// ================ Interface ================
template <typename T, typename PatternType>
void SparsityPattern<T, PatternType>::computePattern() {
    if (pattern != nullptr) return; // Pattern already computed

    if constexpr (std::is_same_v<PatternType, SimplePattern>)
        computeSimplePattern();
    else if constexpr (std::is_same_v<PatternType, GlobalThresholdPattern>)
        computeGlobalThresholdPattern(type.globalThreshold);
    else if constexpr (std::is_same_v<PatternType, ColumnThresholdPattern>)
        computeColumnThresholdPattern(type.columnThreshold);
    else if constexpr (std::is_same_v<PatternType, FixedNNZPattern>)
        computeFixedNNZPattern(type.fixedNNZ);
    else if constexpr (std::is_same_v<PatternType, combinedThresholdPattern>)
        computeCombinedPattern(type.thresh);
    else
        static_assert("Unsupported pattern type");
}

template <typename T, typename PatternType>
inline const CSRMatrix<int> *SparsityPattern<T, PatternType>::getPattern() const {
    assert(pattern != nullptr && "Compute the sparsity pattern first.");
    return pattern;
}

template <typename T, typename PatternType>
inline size_t SparsityPattern<T, PatternType>::getNNZ() const {
    if (!pattern)
        return 0;
    return pattern->nnz;
}

template <typename T, typename PatternType>
SparsityPattern<T, PatternType>::~SparsityPattern() {
    delete pattern;
}

template <typename T, typename PatternType>
void SparsityPattern<T, PatternType>::computeSimplePattern() {
    std::vector<int> patternValues(originalMatrix->nnz, 1);
    pattern = new CSRMatrix<int>(originalMatrix->row_num,
                                               originalMatrix->col_num,
                                               patternValues,
                                               *originalMatrix->row_pointers,
                                               *originalMatrix->col_indices);
}

template <typename T, typename PatternType>
void SparsityPattern<T, PatternType>::computeGlobalThresholdPattern(const double globalThreshold) {
    const std::shared_ptr<std::vector<T>> diagonal = originalMatrix->diagonal(true);

    // Copying the input values for diaonal scaling
    std::vector<T> scaledValues = *(originalMatrix->vals);
    // Perform the diagonal scaling in place
    diagonalScaling(scaledValues, *diagonal);

    // Create the sparsity pattern
    pattern = new CSRMatrix<int>();
    pattern->row_num = originalMatrix->row_num;
    pattern->col_num = originalMatrix->col_num;
    pattern->row_pointers = new std::vector<size_t>(originalMatrix->row_num + 1);

    // Count number of non-zero elements in each row
    if constexpr (SEQUENTIAL) {
        // keep track of the diagonal elements
        std::vector<bool> hasDiagonal(originalMatrix->row_num, false);

        for (size_t i = 0; i < originalMatrix->row_num; ++i) {
            const size_t start = (*(originalMatrix->row_pointers))[i];
            const size_t end = (*(originalMatrix->row_pointers))[i + 1];
            for (size_t j = start; j < end; ++j) {
                const size_t colIdx = (*(originalMatrix->col_indices))[j];
                if (colIdx == i) {
                    hasDiagonal[i] = true;
                    ++(*(pattern->row_pointers))[i + 1];
                } else if (std::abs(scaledValues[j]) > globalThreshold) {
                    ++(*(pattern->row_pointers))[i + 1];
                }
            }

            if (!hasDiagonal[i]) {
                ++(*(pattern->row_pointers))[i + 1];
            }
        }

        // Calculate prefix sum and get non zero count
        // @todo: include it in the scanRowSize function
        pattern->nnz = pattern->scanRowSize();
        pattern->col_indices = new std::vector<size_t>(pattern->nnz);
        pattern->vals = new std::vector<int>(pattern->nnz, 1);

        for (size_t i = 0; i < originalMatrix->row_num; ++i) {
            const size_t start = (*(originalMatrix->row_pointers))[i];
            const size_t end = (*(originalMatrix->row_pointers))[i + 1];
            auto patternColIdx = pattern->col_indices->begin() + (*(pattern->row_pointers))[i];

            for (size_t j = start; j < end; ++j) {
                const size_t colIdx = (*(originalMatrix->col_indices))[j];

                if (!hasDiagonal[i] && colIdx > i) {
                    *patternColIdx++ = i;
                    hasDiagonal[i] = true;
                }

                if ((std::abs(scaledValues[j]) > globalThreshold) || colIdx == i) {
                    *patternColIdx++ = (*(originalMatrix->col_indices))[j];
                }
            }

            if (!hasDiagonal[i]) {
                *patternColIdx++ = i;
            }
        }        
    } else {
        // keep track of the diagonal elements for each thread
        std::vector<std::vector<bool>> hasDiagonal(num_threads, std::vector<bool>(originalMatrix->row_num, false));
        auto fillRowPointer = [&](size_t start, size_t end, int thread_id) {
            for (size_t i = start; i < end; ++i) {
                const size_t rowStart = (*(originalMatrix->row_pointers))[i];
                const size_t rowEnd = (*(originalMatrix->row_pointers))[i + 1];
                for (size_t j = rowStart; j < rowEnd; ++j) {
                    const size_t colIdx = (*(originalMatrix->col_indices))[j];
                    if (colIdx == i) {
                        hasDiagonal[thread_id][i] = true;
                        ++(*(pattern->row_pointers))[i + 1];
                    } else if (std::abs(scaledValues[j]) > globalThreshold) {
                        ++(*(pattern->row_pointers))[i + 1];
                    }
                }

                if (!hasDiagonal[thread_id][i]) {
                    ++(*(pattern->row_pointers))[i + 1];
                }
            }
        };

        launchThreadsWithID(originalMatrix->row_num, fillRowPointer);

        // Calculate prefix sum and get non zero count
        // @todo: include it in the scanRowSize function
        pattern->nnz = pattern->scanRowSize();
        pattern->col_indices = new std::vector<size_t>(pattern->nnz);
        pattern->vals = new std::vector<int>(pattern->nnz, 1);

        auto computeColumnIndices = [&](size_t start, size_t end, int thread_id) {
            for (size_t i = start; i < end; ++i) {
                const size_t rowStart = (*(originalMatrix->row_pointers))[i];
                const size_t rowEnd = (*(originalMatrix->row_pointers))[i + 1];
                auto patternColIdx = pattern->col_indices->begin() + (*(pattern->row_pointers))[i];

                for (size_t j = rowStart; j < rowEnd; ++j) {
                    const size_t colIdx = (*(originalMatrix->col_indices))[j];

                    if (!hasDiagonal[thread_id][i] && colIdx > i) {
                        *patternColIdx++ = i;
                        hasDiagonal[thread_id][i] = true;
                    }

                    if ((std::abs(scaledValues[j]) > globalThreshold) || colIdx == i) {
                        *patternColIdx++ = (*(originalMatrix->col_indices))[j];
                    }
                }

                if (!hasDiagonal[thread_id][i]) {
                    *patternColIdx++ = i;
                }
            }
        };

        launchThreadsWithID(originalMatrix->row_num, computeColumnIndices);
    }
    extend_pattern(*pattern, this->level);
}

// @todo: two loops can be done in parallel
template <typename T, typename PatternType>
void SparsityPattern<T, PatternType>::diagonalScaling(std::vector<T> &values, const std::vector<T> &diagonal) {
    const std::vector<size_t> *rowPointers = originalMatrix->row_pointers;
    const std::vector<size_t> *colIndices = originalMatrix->col_indices;

    // Pre multiplying and post multiplying - (D^-1/2 * A * D^-1/2)
    // multiply each diagonal element with corresponding row in the matrix
    // [a_i] = d_i * [a_i]
    // multiply each diagonal element with corresponding column in the matrix
    // Multiplying each row elements with their corresponding diagonal element (same idx)
    // [a_i]^T = d_i * [a_i]^T
    // @todo: this loop can be done in parallel
    if constexpr (SEQUENTIAL) {
        for (size_t i = 0; i < originalMatrix->row_num; ++i) {
            const size_t start = (*(rowPointers))[i];
            const size_t end = (*(rowPointers))[i + 1];
            for (size_t j = start; j < end; ++j) {
                size_t idx = (*(colIndices))[j];
                values[j] *= diagonal[i] * diagonal[idx]; // diagonal[idx] is a random memory access
            }
        }
    } else {
        auto f = [&](size_t start, size_t end) {
            for (size_t i = start; i < end; ++i) {
                const size_t rowStart = (*(rowPointers))[i];
                const size_t rowEnd = (*(rowPointers))[i + 1];
                for (size_t j = rowStart; j < rowEnd; ++j) {
                    size_t idx = (*(colIndices))[j];
                    values[j] *= diagonal[i] * diagonal[idx]; // diagonal[idx] is a random memory access
                }
            }
        };

        launchThreads(originalMatrix->row_num, f);
    }
}

template <typename T, typename PatternType>
void SparsityPattern<T, PatternType>::computeColumnThresholdPattern(const double tau) {
    pattern = new CSRMatrix<int>();
    pattern->row_num = originalMatrix->row_num;
    pattern->col_num = originalMatrix->col_num;
    pattern->row_pointers = new std::vector<size_t>(originalMatrix->row_num + 1);
    
    // Each thread will only access their partition of the array
    std::vector<T> thresholds(originalMatrix->row_num, 0); // store the threshold value for each row

    // Count number of non-zero elements in each row
    if constexpr (SEQUENTIAL) {
        // keep track of the diagonal elements
        std::vector<bool> hasDiagonal(originalMatrix->row_num, false);
        for (size_t i = 0; i < originalMatrix->row_num; ++i) {
            const size_t start = (*(originalMatrix->row_pointers))[i];
            const size_t end = (*(originalMatrix->row_pointers))[i + 1];
            const T maxVal = *(std::max_element(originalMatrix->vals->begin() + start,
                                                originalMatrix->vals->begin() + end,
                                                [](const T &a, const T &b)
                                                { return std::abs(a) < std::abs(b); }));
            thresholds[i] = (1 - tau) * maxVal;

            for (size_t j = start; j < end; ++j) {
                const size_t colIdx = (*(originalMatrix->col_indices))[j];

                if (colIdx == i) {
                    hasDiagonal[i] = true;
                    ++(*(pattern->row_pointers))[i + 1];
                } else if (std::abs((*(originalMatrix->vals))[j]) > thresholds[i]) {
                    ++(*(pattern->row_pointers))[i + 1];
                }
            }

            if (!hasDiagonal[i]) {
                ++(*(pattern->row_pointers))[i + 1];
            }
        }

        pattern->nnz = pattern->scanRowSize();
        pattern->col_indices = new std::vector<size_t>(pattern->nnz);
        pattern->vals = new std::vector<int>(pattern->nnz, 1);

        for (size_t i = 0; i < originalMatrix->row_num; ++i) {
            const size_t start = (*(originalMatrix->row_pointers))[i];
            const size_t end = (*(originalMatrix->row_pointers))[i + 1];
            size_t patternColIdx = (*(pattern->row_pointers))[i];

            for (size_t j = start; j < end; ++j) {
                const size_t colIdx = (*(originalMatrix->col_indices))[j];

                if (!hasDiagonal[i] && colIdx > i) {
                    (*(pattern->col_indices))[patternColIdx] = i;
                    ++patternColIdx;
                    hasDiagonal[i] = true;
                }

                if (std::abs((*(originalMatrix->vals))[j]) > thresholds[i] || colIdx == i) {
                    (*(pattern->col_indices))[patternColIdx] = (*(originalMatrix->col_indices))[j];
                    ++patternColIdx;
                }
            }

            if (!hasDiagonal[i]) {
                (*(pattern->col_indices))[patternColIdx] = i;
                ++patternColIdx;
            }
        }
    } else {
        // keep track of the diagonal elements
        std::vector<std::vector<bool>> hasDiagonal(num_threads, std::vector<bool>(originalMatrix->row_num, false));
        auto computeRowPointers = [&](size_t start, size_t end, int thread_id) {
            for (size_t i = start; i < end; ++i) {
                const size_t start = (*(originalMatrix->row_pointers))[i];
                const size_t end = (*(originalMatrix->row_pointers))[i + 1];
                const T maxVal = *(std::max_element(originalMatrix->vals->begin() + start,
                                                    originalMatrix->vals->begin() + end,
                                                    [](const T &a, const T &b)
                                                    { return std::abs(a) < std::abs(b); }));
                thresholds[i] = (1 - tau) * maxVal;

                for (size_t j = start; j < end; ++j) {
                    const size_t colIdx = (*(originalMatrix->col_indices))[j];

                    if (colIdx == i) {
                        hasDiagonal[thread_id][i] = true;
                        ++(*(pattern->row_pointers))[i + 1];
                    } else if (std::abs((*(originalMatrix->vals))[j]) > thresholds[i]) {
                        ++(*(pattern->row_pointers))[i + 1];
                    }
                }

                if (!hasDiagonal[thread_id][i]) {
                    ++(*(pattern->row_pointers))[i + 1];
                }
            }
        };

        launchThreadsWithID(originalMatrix->row_num, computeRowPointers);

        pattern->nnz = pattern->scanRowSize();
        pattern->col_indices = new std::vector<size_t>(pattern->nnz);
        pattern->vals = new std::vector<int>(pattern->nnz, 1);

        auto computeColumnIndices = [&](size_t start, size_t end, int thread_id) {
            for (size_t i = start; i < end; ++i) {
                const size_t start = (*(originalMatrix->row_pointers))[i];
                const size_t end = (*(originalMatrix->row_pointers))[i + 1];
                size_t patternColIdx = (*(pattern->row_pointers))[i];

                for (size_t j = start; j < end; ++j) {
                    const size_t colIdx = (*(originalMatrix->col_indices))[j];

                    if (!hasDiagonal[thread_id][i] && colIdx > i) {
                        (*(pattern->col_indices))[patternColIdx] = i;
                        ++patternColIdx;
                        hasDiagonal[thread_id][i] = true;
                    }

                    if (std::abs((*(originalMatrix->vals))[j]) > thresholds[i] || colIdx == i) {
                        (*(pattern->col_indices))[patternColIdx] = (*(originalMatrix->col_indices))[j];
                        ++patternColIdx;
                    }
                }

                if (!hasDiagonal[thread_id][i]) {
                    (*(pattern->col_indices))[patternColIdx] = i;
                    ++patternColIdx;
                }
            }
        };

        launchThreadsWithID(originalMatrix->row_num, computeColumnIndices);
    }

    extend_pattern(*pattern, this->level);
}

template <typename T, typename PatternType>
void SparsityPattern<T, PatternType>::computeFixedNNZPattern(const size_t lfil) {
    pattern = new CSRMatrix<int>();
    pattern->row_num = originalMatrix->row_num;
    pattern->col_num = originalMatrix->col_num;
    pattern->row_pointers = new std::vector<size_t>(originalMatrix->row_num + 1);

    // Array of priority queues, one for each row
    // Each thread will access only its part
    std::vector<std::priority_queue<std::pair<T, size_t>, std::vector<std::pair<T, size_t>>, std::greater<std::pair<T, size_t>>>> lfilLargestElements(pattern->row_num);

    // Count number of non-zero elements in each row
    if constexpr (SEQUENTIAL) {
        for (size_t i = 0; i < originalMatrix->row_num; ++i) {
            const size_t rowStart = (*(originalMatrix->row_pointers))[i];
            const size_t rowEnd = (*(originalMatrix->row_pointers))[i + 1];
            for (size_t j = rowStart; j < rowEnd; ++j) {
                if (lfilLargestElements[i].size() < lfil) {
                    lfilLargestElements[i].push(std::pair<T, size_t>(std::abs((*(originalMatrix->vals))[j]), (*(originalMatrix->col_indices))[j]));
                    ++(*(pattern->row_pointers))[i + 1];
                } else if (lfilLargestElements[i].top().first < std::abs((*(originalMatrix->vals))[j])) {
                    lfilLargestElements[i].pop();
                    lfilLargestElements[i].push(std::pair<T, size_t>(std::abs((*(originalMatrix->vals))[j]), (*(originalMatrix->col_indices))[j]));
                }
            }
        }
    } else {
        auto computeRowPointers = [&](size_t start, size_t end) {
            for (size_t i = start; i < end; ++i) {
                const size_t rowStart = (*(originalMatrix->row_pointers))[i];
                const size_t rowEnd = (*(originalMatrix->row_pointers))[i + 1];
                for (size_t j = rowStart; j < rowEnd; ++j) {
                    if (lfilLargestElements[i].size() < lfil) {
                        lfilLargestElements[i].push(std::pair<T, size_t>(std::abs((*(originalMatrix->vals))[j]), (*(originalMatrix->col_indices))[j]));
                        ++(*(pattern->row_pointers))[i + 1];
                    } else if (lfilLargestElements[i].top().first < std::abs((*(originalMatrix->vals))[j])) {
                        lfilLargestElements[i].pop();
                        lfilLargestElements[i].push(std::pair<T, size_t>(std::abs((*(originalMatrix->vals))[j]), (*(originalMatrix->col_indices))[j]));
                    }
                }
            }
        };

        launchThreads(originalMatrix->row_num, computeRowPointers);
    }

    pattern->nnz = pattern->scanRowSize();
    pattern->col_indices = new std::vector<size_t>(pattern->nnz);
    pattern->vals = new std::vector<int>(pattern->nnz, 1);

    // Compute the column indices
    if constexpr (SEQUENTIAL) {
        for (size_t i = 0; i < originalMatrix->row_num; ++i) {
            const size_t rowStart = (*(pattern->row_pointers))[i];
            const size_t rowEnd = (*(pattern->row_pointers))[i + 1];
            size_t patternColIdx = (*(pattern->row_pointers))[i];

            for (size_t j = rowStart; j < rowEnd; ++j) {
                (*(pattern->col_indices))[patternColIdx] = lfilLargestElements[i].top().second;
                lfilLargestElements[i].pop();
                ++patternColIdx;
            }
        }
    } else {
        auto computeColumnIndices = [&](size_t start, size_t end) {
            for (size_t i = start; i < end; ++i) {
                const size_t rowStart = (*(pattern->row_pointers))[i];
                const size_t rowEnd = (*(pattern->row_pointers))[i + 1];
                size_t patternColIdx = (*(pattern->row_pointers))[i];

                for (size_t j = rowStart; j < rowEnd; ++j) {
                    (*(pattern->col_indices))[patternColIdx] = lfilLargestElements[i].top().second;
                    lfilLargestElements[i].pop();
                    ++patternColIdx;
                }
            }
        };

        launchThreads(originalMatrix->row_num, computeColumnIndices);
    }

    // Sort the column indices in each row
    pattern->sortRows();
    extend_pattern(*pattern, this->level);
}

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