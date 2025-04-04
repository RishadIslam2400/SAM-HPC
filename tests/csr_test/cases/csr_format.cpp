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

    CSRMatrixMock<int> m1(3, 4);
    m1.set(1, 0, 0)
        .set(4, 0, 2)
        .set(5, 0, 3)
        .set(2, 1, 0)
        .set(-1, 1, 1)
        .set(3, 2, 2)
        .set(2, 2, 3);

    std::vector<int> rowPointers1 = {0, 3, 5, 7};
    assertEquals<std::vector<int>>(rowPointers1, m1.getRowPointers(), "Incorrect internal row pointers");

    std::vector<int> colIndices1 = {0, 2, 3, 0, 1, 2, 3};
    assertEquals<std::vector<int>>(colIndices1, m1.getColIndices(), "Incorrect internal column indices");

    std::vector<int> vals1 = {1, 4, 5, 2, -1, 3, 2};
    assertEquals<std::vector<int>>(vals1, m1.getValues(), "Incorrect internal values");

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

    CSRMatrixMock<int> m2(3, 4);
    m2.set(10, 0, 0)
        .set(2, 0, 3)
        .set(3, 2, 0)
        .set(1, 2, 1)
        .set(4, 2, 3);
    
    std::vector<int> rowPointers2 = {0, 2, 2, 5};
    assertEquals<std::vector<int>>(rowPointers2, m2.getRowPointers(), "Incorrect internal row pointers");

    std::vector<int> colIndices2 = {0, 3, 0, 1, 3};
    assertEquals<std::vector<int>>(colIndices2, m2.getColIndices(), "Incorrect internal column indices");

    std::vector<int> vals2 = {10, 2, 3, 1, 4};
    assertEquals<std::vector<int>>(vals2, m2.getValues(), "Incorrect internal values");

    /*
        Previous matrix after adding non-zero element to empty row

        should be stored as
        rows:    [ 0, 2, 3, 6 ]
        columns: [ 0, 3, 1, 0, 1, 3 ]
        values:  [ 10, 2, 5, 3, 1, 4 ]
     */

    CSRMatrixMock<int> m3 = m2;
    m3.set(5, 1, 1);

    std::vector<int> rowPointers3 = {0, 2, 3, 6};
    assertEquals<std::vector<int>>(rowPointers3, m3.getRowPointers(), "Incorrect internal row pointers");

    std::vector<int> colIndices3 = {0, 3, 1, 0, 1, 3};
    assertEquals<std::vector<int>>(colIndices3, m3.getColIndices(), "Incorrect internal column indices");

    std::vector<int> vals3 = {10, 2, 5, 3, 1, 4};
    assertEquals<std::vector<int>>(vals3, m3.getValues(), "Incorrect internal values");

    /*
        Previous matrix with removed the only non-zero element on 2nd row (should be equal to 2nd matrix)

        should be stored as
        rows:    [ 0, 2, 2, 5 ]
        columns: [ 0, 3, 0, 1, 3 ]
        values:  [ 10, 2, 3, 1, 4 ]
     */

    CSRMatrixMock<int> m4 = m3;
    m4.set(0, 1, 1);

    assertEquals<std::vector<int>>(rowPointers2, m4.getRowPointers(), "Incorrect internal row pointers");
    assertEquals<std::vector<int>>(colIndices2, m4.getColIndices(), "Incorrect internal column indices");
    assertEquals<std::vector<int>>(vals2, m4.getValues(), "Incorrect internal values");

    // Create the first matrix using vectors
    CSRMatrixMock<int> m5(3, 4, vals1, rowPointers1, colIndices1);

    assertEquals<std::vector<int>>(rowPointers1, m5.getRowPointers(), "Incorrect internal row pointers");
    assertEquals<std::vector<int>>(colIndices1, m5.getColIndices(), "Incorrect internal column indices");
    assertEquals<std::vector<int>>(vals1, m5.getValues(), "Incorrect internal values");

    std::cout << " OK" << std::endl;
}