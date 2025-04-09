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
    SparseMatrix::CSRMatrix<int> m1(3, 4, vals1, rowPointers1, colIndices1);

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
   SparseMatrix::CSRMatrix<int> m1(3, 4, vals1, rowPointers1, colIndices1);

   // Transpose
   std::shared_ptr<SparseMatrix::CSRMatrix<int>> m1_transpose = m1.transpose();

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
   SparseMatrix::CSRMatrix<int> m2(4, 3, vals2, rowPointers2, colIndices2);
  
   assertEquals<SparseMatrix::CSRMatrix<int>>(m2, *m1_transpose, "Incorrect transpose");

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

       should be stored as
       rows:    [ 0, 3, 5, 7 ]
       columns: [ 0, 2, 3, 0, 1, 2, 3 ]
       values:  [ 1, 4, 5, 2, -1, 4, 2 ]
   */

    // Generate the matrix
    std::vector<size_t> rowPointers1 = {0, 3, 5, 7};
    std::vector<size_t> colIndices1 = {0, 2, 3, 0, 1, 2, 3};
    std::vector<double> vals1 = {1, 4, 5, 2, -1, 4, 2};
    SparseMatrix::CSRMatrix<double> m1(3, 4, vals1, rowPointers1, colIndices1);

    // Diagonal
    std::shared_ptr<std::vector<double>> m1_diag = m1.diagonal();
    std::shared_ptr<std::vector<double>> m1_diag_inverse = m1.diagonal(DiagonalOperation::Inv);
    std::shared_ptr<std::vector<double>> m1_diag_abs = m1.diagonal(DiagonalOperation::Abs);
    std::shared_ptr<std::vector<double>> m1_diag_abs_sqrt = m1.diagonal(DiagonalOperation::AbsSqrt);
    std::shared_ptr<std::vector<double>> m1_diag_abs_inverse = m1.diagonal(DiagonalOperation::AbsInv);
    std::shared_ptr<std::vector<double>> m1_diag_abs_inverse_sqrt = m1.diagonal(DiagonalOperation::AbsInvSqrt);

    std::vector<double> originalDiag = {1, -1, 4};
    std::vector<double> inverseDiag = {1, -1, 0.25};
    std::vector<double> absDiag = {1, 1, 4};
    std::vector<double> absSqrtDiag = {1, 1, 2};
    std::vector<double> absInverseDiag = {1, 1, 0.25};
    std::vector<double> absInverseSqrtDiag = {1, 1, 0.5};

    assertEquals<std::vector<double>>(originalDiag, *m1_diag, "Incorrect diagonal");
    assertEquals<std::vector<double>>(inverseDiag, *m1_diag_inverse, "Incorrect inverse diagonal");
    assertEquals<std::vector<double>>(absDiag, *m1_diag_abs, "Incorrect absolute diagonal");
    assertEquals<std::vector<double>>(absSqrtDiag, *m1_diag_abs_sqrt, "Incorrect absolute sqrt diagonal");
    assertEquals<std::vector<double>>(absInverseDiag, *m1_diag_abs_inverse, "Incorrect absolute inverse diagonal");
    assertEquals<std::vector<double>>(absInverseSqrtDiag, *m1_diag_abs_inverse_sqrt, "Incorrect absolute inverse sqrt diagonal");

    std::cout << "OK" << std::endl;
}