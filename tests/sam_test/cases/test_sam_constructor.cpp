#include "sam.hpp"
#include "testlib.hpp"
#include "helpers.hpp"

#include <iostream>

void testSamConstructor1()
{
    std::cout << "SAM constructor 1..." << std::flush;

    // Create source and target matrices
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
    SparseMatrix::CSRMatrix<int> source(3, 4, vals1, rowPointers1, colIndices1);
    SparseMatrix::CSRMatrix<int> target(source);
    SparsityPattern<int, GlobalThresholdPattern> pattern(source, GlobalThresholdPattern{0.001});

    // Create the SAM object
    SparseApproximateMap<int, GlobalThresholdPattern> map(target, source, pattern);
    
    /*
    assertEquals<SparseMatrix::CSRMatrix<int>>(*map.targetMatrix, target);
    assertEquals<SparseMatrix::CSRMatrix<int>>(*map.sourceMatrix, source);
    assertEquals<SparseMatrix::CSRMatrix<int>*>(map.mappingMatrix, nullptr);
    assertEquals<SparsityPattern<int, SimplePattern>>(*map.sparsityPattern, pattern);
    */

    std::cout << "OK" << std::endl;
}

void testSamConstructor2()
{
    std::cout << "SAM constructor 2..." << std::flush;

    // Create source and target matrices
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
    SparseMatrix::CSRMatrix<int> source(3, 4, vals1, rowPointers1, colIndices1);
    SparseMatrix::CSRMatrix<int> target(source);
    GlobalThresholdPattern type{0.001};
    SparsityPattern<int, GlobalThresholdPattern> pattern(source, type);

    // Create the SAM object
    SparseApproximateMap<int, GlobalThresholdPattern> map(target, source, type);

    /*
    assertEquals<SparseMatrix::CSRMatrix<int>>(*map.targetMatrix, target);
    assertEquals<SparseMatrix::CSRMatrix<int>>(*map.sourceMatrix, source);
    assertEquals<SparsityPattern<int, GlobalThresholdPattern>>(*map.sparsityPattern, pattern);
    assertEquals<SparseMatrix::CSRMatrix<int>*>(map.mappingMatrix, nullptr);
    */

    std::cout << "OK" << std::endl;
}