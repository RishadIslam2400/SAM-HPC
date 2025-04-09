#include "CSRMatrix.hpp"

enum class SparsityPatternType
{
    SIMPLE,
    GLOBAL_THRESH,
    COLUMN_THRESH,
    FIXED_NNZ
};

struct SparsityPatternParams
{
    double globalThreshold;
    double columnThreshold;
    size_t fixedNNZ;
};

template <typename T>
class SparsityPattern
{
private:
    SparseMatrix::CSRMatrix<int> *pattern;
    const SparseMatrix::CSRMatrix<T> *originalMatrix;
    SparsityPatternType type;

    // ================ Sparsity Pattern Computation ================
    void computeSimplePattern();
    void computeGlobalThresholdPattern(const double globalThreshold);
    void computeColumnThresholdPattern();
    void computeFixedNNZPattern();

    // =============== Helper Functions ================
    void diagonalScaling(std::vector<T> &values, const std::vector<T> &diagonal);

public:
    SparsityPattern() = delete;
    SparsityPattern(SparseMatrix::CSRMatrix<T> &originalMatrix, SparsityPatternType type);
    void computePattern(const SparsityPatternParams &params = {});

    const SparseMatrix::CSRMatrix<int> *getPattern() const;

    ~SparsityPattern();
};

template <typename T>
void SparsityPattern<T>::computeSimplePattern()
{
    std::vector<int> patternValues(originalMatrix->nnz, 1);
    pattern = new SparseMatrix::CSRMatrix<int>(originalMatrix->row_num,
                                               originalMatrix->col_num,
                                               patternValues,
                                               *originalMatrix->row_pointers,
                                               *originalMatrix->col_indices);
}

template <typename T>
void SparsityPattern<T>::computeGlobalThresholdPattern(const double globalThreshold)
{
    const std::shared_ptr<std::vector<T>> diagonal = originalMatrix->diagonal(DiagonalOperation::AbsInvSqrt);

    // Copying the input values for diaonal scaling
    std::vector<T> scaledValues = *(originalMatrix->vals);
    // Perform the diagonal scaling in place
    diagonalScaling(scaledValues, *diagonal);

    // Create the sparsity pattern
    pattern = new SparseMatrix::CSRMatrix<int>();
    pattern->row_num = originalMatrix->row_num;
    pattern->col_num = originalMatrix->col_num;
    pattern->row_pointers = new std::vector<size_t>(originalMatrix->row_num + 1);

    // Count number of non-zero elements in each row
    for (size_t i = 0; i < originalMatrix->row_num; ++i)
    {
        const size_t start = (*(originalMatrix->row_pointers))[i];
        const size_t end = (*(originalMatrix->row_pointers))[i + 1];
        for (size_t j = start; j < end; ++j)
        {
            if (std::abs(scaledValues[j]) > globalThreshold)
            {
                ++(*(pattern->row_pointers))[i + 1];
            }
        }
    }

    // Calculate prefix sum and get non zero count
    size_t nonZeroCount = pattern->scanRowSize();
    pattern->col_indices = new std::vector<size_t>(nonZeroCount);
    pattern->vals = new std::vector<int>(nonZeroCount, 1);

    // Compute the column indices
    for (size_t i = 0; i < originalMatrix->row_num; ++i)
    {
        const size_t start = (*(originalMatrix->row_pointers))[i];
        const size_t end = (*(originalMatrix->row_pointers))[i + 1];
        size_t patternColIdx = (*(pattern->row_pointers))[i];

        for (size_t j = start; j < end; ++j)
        {
            if ((std::abs(scaledValues[j]) > globalThreshold))
            {
                (*(pattern->col_indices))[patternColIdx] = (*(originalMatrix->col_indices))[j];
                ++patternColIdx;
            }
        }
    }
}

// @todo: two loops can be done in parallel
template <typename T>
void SparsityPattern<T>::diagonalScaling(std::vector<T> &values, const std::vector<T> &diagonal)
{
    const std::vector<size_t> *rowPointers = originalMatrix->row_pointers;
    const std::vector<size_t> *colIndices = originalMatrix->col_indices;

    // Pre multiplying and post multiplying - (D^-1/2 * A * D^-1/2)
    // multiply each diagonal element with corresponding row in the matrix
    // [a_i] = d_i * [a_i]
    // multiply each diagonal element with corresponding column in the matrix
    // Multiplying each row elements with their corresponding diagonal element (same idx)
    // [a_i]^T = d_i * [a_i]^T
    // @todo: this loop can be done in parallel
    for (size_t i = 0; i < originalMatrix->row_num; ++i)
    {
        const size_t start = (*(rowPointers))[i];
        const size_t end = (*(rowPointers))[i + 1];
        for (size_t j = start; j < end; ++j)
        {
            size_t idx = (*(colIndices))[j];
            values[j] *= diagonal[i] * diagonal[idx]; // diagonal[idx] is a random memory access
        }
    }
}

template <typename T>
inline void SparsityPattern<T>::computeColumnThresholdPattern()
{
}

template <typename T>
inline void SparsityPattern<T>::computeFixedNNZPattern()
{
}

template <typename T>
SparsityPattern<T>::SparsityPattern(SparseMatrix::CSRMatrix<T> &originalMatrix, SparsityPatternType t)
    : originalMatrix(&originalMatrix), type(t)
{
    pattern = nullptr;
}

template <typename T>
void SparsityPattern<T>::computePattern(const SparsityPatternParams &params)
{
    assert(pattern == nullptr && "Sparsity pattern already computed.");

    switch (type)
    {
    case SparsityPatternType::SIMPLE:
        computeSimplePattern();
        break;
    case SparsityPatternType::GLOBAL_THRESH:
        computeGlobalThresholdPattern(params.globalThreshold);
        break;
    case SparsityPatternType::COLUMN_THRESH:
        computeColumnThresholdPattern();
        break;
    case SparsityPatternType::FIXED_NNZ:
        computeFixedNNZPattern();
        break;
    default:
        break;
    }
}

template <typename T>
inline const SparseMatrix::CSRMatrix<int> *SparsityPattern<T>::getPattern() const
{
    assert(pattern != nullptr && "Sparsity pattern not computed.");
    return pattern;
}

template <typename T>
SparsityPattern<T>::~SparsityPattern()
{
    delete pattern;
}