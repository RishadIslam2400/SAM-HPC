#include "CSRMatrix.hpp"
#include "testlib.hpp"

#include <memory>

void testSortRows()
{
    std::cout << "sort rows..." << std::flush;
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

    // rearrange the columns and values
    std::vector<size_t> rowPointers1 = {0, 3, 5, 7};
    std::vector<size_t> colIndices1 = {2, 3, 0, 1, 0, 3, 2};
    std::vector<int> vals1 = {4, 5, 1, -1, 2, 2 , 3};
    CSRMatrix<int> m1(3, 4, vals1, rowPointers1, colIndices1);

    std::vector<size_t> rowPointersCorrect = {0, 3, 5, 7};
    std::vector<size_t> colIndicesCorrect = {0, 2, 3, 0, 1, 2, 3};
    std::vector<int> valsCorrect = {1, 4, 5, 2, -1, 3, 2};

    m1.sortRows();
    assertEquals<std::vector<size_t>>(rowPointersCorrect, *(m1.row_pointers), "Incorrect internal row pointers");
    assertEquals<std::vector<size_t>>(colIndicesCorrect, *(m1.col_indices), "Incorrect internal column indices");
    assertEquals<std::vector<int>>(valsCorrect, *(m1.vals), "Incorrect internal values");

    std::cout << "OK" << std::endl;
}

void testTranspose()
{
    std::cout << "transpose..." << std::flush;
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
   CSRMatrix<int> m1(3, 4, vals1, rowPointers1, colIndices1);

   // Transpose
   std::shared_ptr<CSRMatrix<int>> m1_transpose = m1.transpose();

   /*
       "Transpose" matrix
       [ 1  2 0 ]
       [ 0 -1 0 ]
       [ 4  0 3 ]
       [ 5  0 2 ]

       should be stored as
       rows:    [ 0, 2, 3, 5, 7 ]
       columns: [ 0, 1, 1, 0, 2, 0, 2 ]
       values:  [ 1, 2, -1, 4, 3, 5, 2 ]
   */

   std::vector<size_t> rowPointers2 = {0, 2, 3, 5, 7};
   std::vector<size_t> colIndices2 = {0, 1, 1, 0, 2, 0, 2};
   std::vector<int> vals2 = {1, 2, -1, 4, 3, 5, 2};
   CSRMatrix<int> m2(4, 3, vals2, rowPointers2, colIndices2);
  
   assertEquals<CSRMatrix<int>>(m2, *m1_transpose, "Incorrect transpose");

   std::cout << "OK" << std::endl;
}

void testDiagonal()
{
    std::cout << "diagonal..." << std::flush;

    /*
       "Standard" matrix
       [ 1  0 4 5 ]
       [ 2 -1 0 0 ]
       [ 0  0 4 2 ]
       [ 0  0 0 0 ]

       should be stored as
       rows:    [ 0, 3, 5, 7, 7 ]
       columns: [ 0, 2, 3, 0, 1, 2, 3 ]
       values:  [ 1, 4, 5, 2, -1, 4, 2 ]
   */

    // Generate the matrix
    std::vector<size_t> rowPointers1 = {0, 3, 5, 7, 7};
    std::vector<size_t> colIndices1 = {0, 2, 3, 0, 1, 2, 3};
    std::vector<double> vals1 = {1, 4, 5, 2, -1, 4, 2};
    CSRMatrix<double> m1(4, 4, vals1, rowPointers1, colIndices1);

    // Diagonal
    std::shared_ptr<std::vector<double>> m1_diag = m1.diagonal();
    std::shared_ptr<std::vector<double>> m1_diag_without_scaling = m1.diagonal(false);
    std::shared_ptr<std::vector<double>> m1_diag_for_scaling = m1.diagonal(true);

    std::vector<double> originalDiag = {1, -1, 4, 0};
    std::vector<double> forScalingDiag = {1, 1, 0.5, 1};

    assertEquals<std::vector<double>>(*m1_diag, originalDiag, "Incorrect diagonal");
    assertEquals<std::vector<double>>(*m1_diag_without_scaling, originalDiag, "Incorrect diagonal");
    assertEquals<std::vector<double>>(*m1_diag_for_scaling, forScalingDiag, "Incorrect diagonal for scaling");

    std::cout << "OK" << std::endl;
}