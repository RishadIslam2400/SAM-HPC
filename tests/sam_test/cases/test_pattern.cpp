#include "CSRMatrix.hpp"
#include "sparsityPattern.hpp"
#include "testlib.hpp"
#include "helpers.hpp"

#include <iostream>

void testSimpleSparsityPattern()
{
    std::cout << "Simple sparsity pattern..." << std::flush;

    /*
        "Standard" matrix
        [ 1  0 4 5 ]
        [ 2 -1 0 0 ]
        [ 0  0 3 2 ]

        should be stored as
        rows:    [ 0, 3, 5, 7 ]
        columns: [ 0, 2, 3, 0, 1, 2, 3 ]
        values:  [ 1, 4, 5, 2, -1, 3, 2 ]
    */

    // Generate the matrix
    std::vector<size_t> rowPointers1 = {0, 3, 5, 7};
    std::vector<size_t> colIndices1 = {0, 2, 3, 0, 1, 2, 3};
    std::vector<int> vals1 = {1, 4, 5, 2, -1, 3, 2};
    SparseMatrix::CSRMatrix<int> m1(3, 4, vals1, rowPointers1, colIndices1);

    // Generate the sparsity pattern
    SparsityPattern<int> simplePattern(m1, SparsityPatternType::SIMPLE);
    simplePattern.computePattern();
    const SparseMatrix::CSRMatrix<int>* recievedPattern = simplePattern.getPattern();

    std::vector<int> patternValuesCorrect(7, 1);

    assertEquals<std::vector<size_t>>(*(recievedPattern->row_pointers), rowPointers1, "Incorrect internal row pointers");
    assertEquals<std::vector<size_t>>(*(recievedPattern->col_indices), colIndices1, "Incorrect internal column indices");
    assertEquals<std::vector<int>>(*(recievedPattern->vals), patternValuesCorrect, "Incorrect internal values");

    std::cout << "OK" << std::endl;
}

void testGlobalSparsityPattern()
{
    std::cout << "Global sparsity pattern..." << std::flush;

    /*
        "Standard" matrix
        [ 1  0 4 5 ]
        [ 2 -1 0 0 ]
        [ 0  0 3 2 ]
        [ 7  0 2 0 ]

        should be stored as
        rows:    [ 0, 3, 5, 7, 9 ]
        columns: [ 0, 2, 3, 0, 1, 2, 3, 0, 2 ]
        values:  [ 1, 4, 5, 2, -1, 3, 2, 7, 2 ]
    */

    // Generate the matrix
    std::vector<size_t> rowPointers1 = {0, 3, 5, 7, 9};
    std::vector<size_t> colIndices1 = {0, 2, 3, 0, 1, 2, 3, 0, 2};
    std::vector<double> vals1 = {1, 4, 5, 2, -1, 3, 2, 7, 2};
    SparseMatrix::CSRMatrix<double> m1(4, 4, vals1, rowPointers1, colIndices1);

    // Generate the sparsity pattern
    SparsityPattern<double> simplePattern(m1, SparsityPatternType::GLOBAL_THRESH);
    SparsityPatternParams params;
    params.globalThreshold = 0.001;
    simplePattern.computePattern(params);
    const SparseMatrix::CSRMatrix<int> *recievedPattern = simplePattern.getPattern();

    std::vector<int> patternValuesCorrect = {1, 1, 1, 1, 1};
    std::vector<size_t> patternRowPointersCorrect = {0, 2, 4, 5, 5};
    std::vector<size_t> patternColIndicesCorrect = {0, 2, 0, 1, 2};

    assertEquals<std::vector<size_t>>(*(recievedPattern->row_pointers), patternRowPointersCorrect, "Incorrect internal row pointers");
    assertEquals<std::vector<size_t>>(*(recievedPattern->col_indices), patternColIndicesCorrect, "Incorrect internal column indices");
    assertEquals<std::vector<int>>(*(recievedPattern->vals), patternValuesCorrect, "Incorrect internal values"); 

    std::cout << "OK" << std::endl;
}