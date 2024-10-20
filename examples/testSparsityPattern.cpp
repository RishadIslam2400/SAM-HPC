#include <iostream>
#include <vector>
#include <iomanip>

#include "sparsityPattern.hpp"
#include "samFunc.hpp"

int main() {
    std::vector<double> values = 
    {1.1, 2.2, 3.3, 4.4, 5.5, 6.6, 7.7, 8.8, 9.9, 10.0, 11.1, 12.2, 13.3, 14.4, 15.5, 16.6, 17.7, 18.8};
    std::vector<size_t> rowIndices = {0, 1, 1, 5, 6, 3, 4, 5, 7, 3, 3, 5, 0, 6, 8, 6, 4, 7};
    std::vector<size_t> colPointers = {0, 2, 5, 9, 10, 12, 12, 14, 15, 16, 18};

    size_t numRows = 10;
    size_t numCols = 10;
    size_t nnz = 18;

    csc_matrix A(values, rowIndices, colPointers, numCols, numRows, nnz);
    csc_matrix S;

    sparsity_pattern_global_thresh(A, 0.5, S);
}