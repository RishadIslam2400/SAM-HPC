#include "CSRMatrixMock.hpp"
#include "testlib.hpp"

void _constructorFail1()
{
    // Create a matrix
    /*
            [0 2 0 1]                                    
        A = [1 3 0 0]                                
            [0 7 4 0]                                    
        row pointers:   [0, 2, 4, 6]                     
        column indices: [1, 3, 0, 1, 1, 2]
        values:         [2, 1, 1, 3, 7, 4]    
                                              
    */
    std::vector<size_t> rowPointersA = {0, 2, 4, 6};
    std::vector<size_t> colIndicesA = {1, 3, 0, 1, 1, 2};
    std::vector<int> valsA = {2, 1, 1, 3, 7, 4};
    CSRMatrix<int> A(0, 0, valsA, rowPointersA, colIndicesA); // size zero
}

void testConstructorFail1()
{
    std::cout << "constructor fail #1..." << std::flush;
    assertException("InvalidDimensionsException", _constructorFail1);
    std::cout << "OK" << std::endl;
}

void _constructorFail2()
{
    // Create a matrix
    /*
            [0 2 0 1]
        A = [1 3 0 0]
            [0 7 4 0]
        row pointers:   [0, 2, 4, 6]
        column indices: [1, 3, 0, 1, 1, 2]
        values:         [2, 1, 1, 3, 7, 4]

    */
    std::vector<size_t> rowPointersA = {0, 2, 4, 6, 8};
    std::vector<size_t> colIndicesA = {1, 3, 0, 1, 1, 2};
    std::vector<int> valsA = {2, 1, 1, 3, 7, 4};
    CSRMatrix<int> A(3, 4, valsA, rowPointersA, colIndicesA); // row pointers out of bounds
}

void testConstructorFail2()
{
    std::cout << "constructor fail #2..." << std::flush;
    assertException("InvalidDimensionsException", _constructorFail2);
    std::cout << " OK" << std::endl;
}

void _constructorFail3()
{
    // Create a matrix
    /*
            [0 2 0 1]
        A = [1 3 0 0]
            [0 7 4 0]
        row pointers:   [0, 2, 4, 6]
        column indices: [1, 3, 0, 1, 1, 2]
        values:         [2, 1, 1, 3, 7, 4]

    */
    std::vector<size_t> rowPointersA = {0, 2, 4, 6};
    std::vector<size_t> colIndicesA = {1, 3, 0, 1, 1, 2, 1};
    std::vector<int> valsA = {2, 1, 1, 3, 7, 4};
    CSRMatrix<int> A(3, 4, valsA, rowPointersA, colIndicesA); // column indices out of bounds
}

void testConstructorFail3()
{
    std::cout << "constructor fail #3..." << std::flush;
    assertException("InvalidDimensionsException", _constructorFail3);
    std::cout << " OK" << std::endl;
}

void _constructorFail4()
{
    // Create a matrix
    /*
            [0 2 0 1]
        A = [1 3 0 0]
            [0 7 4 0]
        row pointers:   [0, 2, 4, 6]
        column indices: [1, 3, 0, 1, 1, 2]
        values:         [2, 1, 1, 3, 7, 4]

    */
    std::vector<size_t> rowPointersA = {0, 2, 4, 6};
    std::vector<size_t> colIndicesA = {1, 3, 0, 1, 1, 2};
    std::vector<int> valsA = {2, 1, 1, 3, 7, 4, 10};
    CSRMatrix<int> A(3, 4, valsA, rowPointersA, colIndicesA); // values out of bounds
}

void testConstructorFail4()
{
    std::cout << "constructor fail #4..." << std::flush;
    assertException("InvalidDimensionsException", _constructorFail4);
    std::cout << " OK" << std::endl;
}

void _constructorFail5()
{
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
    std::vector<size_t> rowPointers1 = {0, 3, 5, 10};
    std::vector<size_t> colIndices1 = {0, 2, 3, 0, 1, 2, 3};
    std::vector<int> vals1 = {1, 4, 5, 2, -1, 3, 2};
    CSRMatrix<int> m1(3, 4, vals1, rowPointers1, colIndices1); // non zero count does not match
}

void testConstructorFail5()
{
    std::cout << "constructor fail #5..." << std::flush;
    assertException("InvalidDimensionsException", _constructorFail5);
    std::cout << " OK" << std::endl;
}

void _constructorFail6()
{
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
    std::vector<size_t> rowPointers1 = {0, 3, 5, 7};
    std::vector<size_t> colIndices1 = {0, 2, 3, 0, 1, 2, 3};
    std::vector<int> vals1 = {1, 4, 5, 2, -1, 3, 2};
    CSRMatrix<int> m1(3, 4, 10, vals1, rowPointers1, colIndices1); // nnz does not match
}

void testConstructorFail6()
{
    std::cout << "constructor fail #6..." << std::flush;
    assertException("InvalidDimensionsException", _constructorFail6);
    std::cout << " OK" << std::endl;
}

void _constructorFail7()
{
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
    std::vector<size_t> rowPointers1 = {0, 3, 5, 7};
    std::vector<size_t> colIndices1 = {0, 2, 3, 0, 1, 2, 3, 8};
    std::vector<int> vals1 = {1, 4, 5, 2, -1, 3};
    CSRMatrix<int> m1(vals1, rowPointers1, colIndices1); // column indices out of bounds
}

void testConstructorFail7()
{
    std::cout << "constructor fail #7..." << std::flush;
    assertException("InvalidDimensionsException", _constructorFail7);
    std::cout << " OK" << std::endl;
}

void _constructorFail8()
{
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
    std::vector<size_t> rowPointers1 = {0, 3, 5, 7};
    std::vector<size_t> colIndices1 = {0, 2, 3, 0, 1, 2, 3};
    std::vector<int> vals1 = {1, 4, 5, 2, -1, 3};
    CSRMatrix<int> m1(vals1, rowPointers1, colIndices1); // values does not match nnz count
}

void testConstructorFail8()
{
    std::cout << "constructor fail #8..." << std::flush;
    assertException("InvalidDimensionsException", _constructorFail8);
    std::cout << " OK" << std::endl;
}

void _constructorFail9()
{
    size_t* rowPointers1 = nullptr;
    size_t* colIndices1 = nullptr;
    int* vals1 = nullptr;
    CSRMatrix<int> m1(3, 4, vals1, rowPointers1, colIndices1); // passing null pointers
}

void testConstructorFail9()
{
    std::cout << "constructor fail #9..." << std::flush;
    assertException("InvalidDimensionsException", _constructorFail9);
    std::cout << " OK" << std::endl;
}