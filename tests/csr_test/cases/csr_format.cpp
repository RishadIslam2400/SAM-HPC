#include "CSRMatrixMock.hpp"
#include "testlib.hpp"
#include "helpers.hpp"

void testInternalStorage()
{
    std::cout << "internal storage..." << std::flush;

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
    CSRMatrixMock<int> m1(3, 4, vals1, rowPointers1, colIndices1);

    assertEquals<std::vector<size_t>>(rowPointers1, *(m1.row_pointers), "Incorrect internal row pointers");
    assertEquals<std::vector<size_t>>(colIndices1, *(m1.col_indices), "Incorrect internal column indices");
    assertEquals<std::vector<int>>(vals1, *(m1.vals), "Incorrect internal values");

    /*
        Matrix with empty row
        [ 10 0 0 2 ]
        [  0 0 0 0 ]
        [  3 1 0 4 ]

        should be stored as
        rows:    [ 0, 2, 2, 5 ]
        columns: [ 0, 3, 0, 1, 3 ]
        values:  [ 10, 2, 3, 1, 4 ]
    */

    std::vector<size_t> rowPointers2 = {0, 2, 2, 5};
    std::vector<size_t> colIndices2 = {0, 3, 0, 1, 3};
    std::vector<int> vals2 = {10, 2, 3, 1, 4};
    CSRMatrixMock<int> m2(3, 4, vals2, rowPointers2, colIndices2);
    
    assertEquals<std::vector<size_t>>(rowPointers2, *(m2.row_pointers), "Incorrect internal row pointers");
    assertEquals<std::vector<size_t>>(colIndices2, *(m2.col_indices), "Incorrect internal column indices");
    assertEquals<std::vector<int>>(vals2, *(m2.vals), "Incorrect internal values");

    /*
        Previous matrix after adding non-zero element to empty row

        should be stored as
        rows:    [ 0, 2, 3, 6 ]
        columns: [ 0, 3, 1, 0, 1, 3 ]
        values:  [ 10, 2, 5, 3, 1, 4 ]
     */

    CSRMatrixMock<int> m3 = m2;
    m3.set(5, 1, 1);

    std::vector<size_t> rowPointers3 = {0, 2, 3, 6};
    assertEquals<std::vector<size_t>>(rowPointers3, *(m3.row_pointers), "Incorrect internal row pointers");

    std::vector<size_t> colIndices3 = {0, 3, 1, 0, 1, 3};
    assertEquals<std::vector<size_t>>(colIndices3, *(m3.col_indices), "Incorrect internal column indices");

    std::vector<int> vals3 = {10, 2, 5, 3, 1, 4};
    assertEquals<std::vector<int>>(vals3, *(m3.vals), "Incorrect internal values");

    /*
        Previous matrix with removed the only non-zero element on 2nd row (should be equal to 2nd matrix)

        should be stored as
        rows:    [ 0, 2, 2, 5 ]
        columns: [ 0, 3, 0, 1, 3 ]
        values:  [ 10, 2, 3, 1, 4 ]
     */

    CSRMatrixMock<int> m4 = m3;
    m4.set(0, 1, 1);

    assertEquals<std::vector<size_t>>(rowPointers2, *(m4.row_pointers), "Incorrect internal row pointers");
    assertEquals<std::vector<size_t>>(colIndices2, *(m4.col_indices), "Incorrect internal column indices");
    assertEquals<std::vector<int>>(vals2, *(m4.vals), "Incorrect internal values");


    // Construct from dense 2D matrix
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

    std::vector<std::vector<int>> denseMatrix = {
        {1, 0, 4, 5},
        {2, -1, 0, 0},
        {0, 0, 3, 2}
    };

    CSRMatrixMock<int> m5(denseMatrix);
    assertEquals<std::vector<size_t>>(rowPointers1, *(m5.row_pointers), "Incorrect internal row pointers");
    assertEquals<std::vector<size_t>>(colIndices1, *(m5.col_indices), "Incorrect internal column indices");
    assertEquals<std::vector<int>>(vals1, *(m5.vals), "Incorrect internal values");

    auto row0 = m5.rowBegin(0);
    assertEquals<size_t>(0, row0.col(), "Incorrect row iterator.");
    assertEquals<int>(1, row0.value(), "Incorrect row iterator.");
    ++row0;
    assertEquals<size_t>(2, row0.col(), "Incorrect row iterator.");
    assertEquals<int>(4, row0.value(), "Incorrect row iterator.");
    ++row0;
    assertEquals<size_t>(3, row0.col(), "Incorrect row iterator.");
    assertEquals<int>(5, row0.value(), "Incorrect row iterator.");

    auto row1 = m5.rowBegin(1);
    assertEquals<size_t>(0, row1.col(), "Incorrect row iterator.");
    assertEquals<int>(2, row1.value(), "Incorrect row iterator.");
    ++row1;
    assertEquals<size_t>(1, row1.col(), "Incorrect row iterator.");
    assertEquals<int>(-1, row1.value(), "Incorrect row iterator.");

    // Construct a matrix using only vectors
    CSRMatrixMock<int> m6(vals1, rowPointers1, colIndices1);
    assertEquals<std::vector<size_t>>(rowPointers1, *(m6.row_pointers), "Incorrect internal row pointers");
    assertEquals<std::vector<size_t>>(colIndices1, *(m6.col_indices), "Incorrect internal column indices");
    assertEquals<std::vector<int>>(vals1, *(m6.vals), "Incorrect internal values");

    // Construct a matrix using non zero count
    CSRMatrixMock<int> m7(3, 4, 7, vals1, rowPointers1, colIndices1);
    assertEquals<std::vector<size_t>>(rowPointers1, *(m7.row_pointers), "Incorrect internal row pointers");
    assertEquals<std::vector<size_t>>(colIndices1, *(m7.col_indices), "Incorrect internal column indices");
    assertEquals<std::vector<int>>(vals1, *(m7.vals), "Incorrect internal values");

    // Construct a matrix using array pointers
    CSRMatrixMock<int> m8(3, 4, vals1.data(), rowPointers1.data(), colIndices1.data());
    assertEquals<std::vector<size_t>>(rowPointers1, *(m8.row_pointers), "Incorrect internal row pointers");
    assertEquals<std::vector<size_t>>(colIndices1, *(m8.col_indices), "Incorrect internal column indices");
    assertEquals<std::vector<int>>(vals1, *(m8.vals), "Incorrect internal values");

    // Constrcut a matrix using array pointers and non zero count
    CSRMatrixMock<int> m9(3, 4, 7, vals1.data(), rowPointers1.data(), colIndices1.data());
    assertEquals<std::vector<size_t>>(rowPointers1, *(m9.row_pointers), "Incorrect internal row pointers");
    assertEquals<std::vector<size_t>>(colIndices1, *(m9.col_indices), "Incorrect internal column indices");
    assertEquals<std::vector<int>>(vals1, *(m9.vals), "Incorrect internal values");

    std::cout << " OK" << std::endl;
}