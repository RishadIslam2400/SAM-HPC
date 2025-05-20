#include "sam.hpp"
#include "sparsityPattern.hpp"
#include "read_mat.hpp"
#include "testlib.hpp"
#include "helpers.hpp"


void testSAMSanityCheck1()
{
    std::cout << "SAM Sanity Check 1..." << std::flush;

    const size_t targetRows = 5;
    const size_t targetCols = 5;
    const size_t targetNNZ = 7;
    const std::vector<size_t> targetRowPointers = {0, 2, 3, 5, 6, 7};
    const std::vector<size_t> targetColIndices = {0, 2, 1, 0, 4, 3, 2};
    const std::vector<double> targetVals = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0};
    const CSRMatrix<double> targetMatrix(targetRows, targetCols, targetNNZ, targetVals, targetRowPointers, targetColIndices);

    const size_t sourceRows = 5;
    const size_t sourceCols = 5;
    const size_t sourceNNZ = 8;
    const std::vector<size_t> sourceRowPointers = {0, 1, 3, 5, 7, 8};
    const std::vector<size_t> sourceColIndices = {1, 0, 4, 1, 2, 0, 3, 2};
    const std::vector<double> sourceVals = {10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0};
    const CSRMatrix<double> sourceMatrix(sourceRows, sourceCols, sourceNNZ, sourceVals, sourceRowPointers, sourceColIndices);

    SparsityPattern<double, SimplePattern> pattern(sourceMatrix, SimplePattern());
    pattern.computePattern();

    SparseApproximateMap<double, SimplePattern> testSAM(targetMatrix, sourceMatrix, pattern);
    testSAM.computeMap();

    std::cout << "OK" << std::endl;
}

void testSAMSanityCheck2()
{
    std::cout << "SAM Sanity Check 2..." << std::flush;

    const size_t targetRows = 5;
    const size_t targetCols = 5;
    const size_t targetNNZ = 7;
    const std::vector<size_t> targetRowPointers = {0, 2, 3, 5, 6, 7};
    const std::vector<size_t> targetColIndices = {0, 2, 1, 0, 4, 3, 2};
    const std::vector<double> targetVals = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0};
    const CSRMatrix<double> targetMatrix(targetRows, targetCols, targetNNZ, targetVals, targetRowPointers, targetColIndices);

    const size_t sourceRows = 5;
    const size_t sourceCols = 5;
    const size_t sourceNNZ = 8;
    const std::vector<size_t> sourceRowPointers = {0, 1, 3, 5, 7, 8};
    const std::vector<size_t> sourceColIndices = {1, 0, 4, 1, 2, 0, 3, 2};
    const std::vector<double> sourceVals = {10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0};
    const CSRMatrix<double> sourceMatrix(sourceRows, sourceCols, sourceNNZ, sourceVals, sourceRowPointers, sourceColIndices);

    GlobalThresholdPattern thresh{0.001};
    SparsityPattern<double, GlobalThresholdPattern> pattern(sourceMatrix, thresh);
    pattern.computePattern();

    SparseApproximateMap<double, GlobalThresholdPattern> testSAM(targetMatrix, sourceMatrix, pattern);
    testSAM.computeMap();

    std::cout << "OK" << std::endl;
}

void testSAMSanityCheck3()
{
    std::cout << "SAM Sanity Check 3..." << std::flush;

    const size_t targetRows = 5;
    const size_t targetCols = 5;
    const size_t targetNNZ = 7;
    const std::vector<size_t> targetRowPointers = {0, 2, 3, 5, 6, 7};
    const std::vector<size_t> targetColIndices = {0, 2, 1, 0, 4, 3, 2};
    const std::vector<double> targetVals = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0};
    const CSRMatrix<double> targetMatrix(targetRows, targetCols, targetNNZ, targetVals, targetRowPointers, targetColIndices);

    const size_t sourceRows = 5;
    const size_t sourceCols = 5;
    const size_t sourceNNZ = 8;
    const std::vector<size_t> sourceRowPointers = {0, 1, 3, 5, 7, 8};
    const std::vector<size_t> sourceColIndices = {1, 0, 4, 1, 2, 0, 3, 2};
    const std::vector<double> sourceVals = {10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0};
    const CSRMatrix<double> sourceMatrix(sourceRows, sourceCols, sourceNNZ, sourceVals, sourceRowPointers, sourceColIndices);

    ColumnThresholdPattern thresh{0.9};
    SparsityPattern<double, ColumnThresholdPattern> pattern(sourceMatrix, thresh);
    pattern.computePattern();

    SparseApproximateMap<double, ColumnThresholdPattern> testSAM(targetMatrix, sourceMatrix, pattern);
    testSAM.computeMap();

    std::cout << "OK" << std::endl;
}

void testSAMSanityCheck4()
{
    std::cout << "SAM Sanity Check 4..." << std::flush;

    const size_t targetRows = 5;
    const size_t targetCols = 5;
    const size_t targetNNZ = 7;
    const std::vector<size_t> targetRowPointers = {0, 2, 3, 5, 6, 7};
    const std::vector<size_t> targetColIndices = {0, 2, 1, 0, 4, 3, 2};
    const std::vector<double> targetVals = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0};
    const CSRMatrix<double> targetMatrix(targetRows, targetCols, targetNNZ, targetVals, targetRowPointers, targetColIndices);

    const size_t sourceRows = 5;
    const size_t sourceCols = 5;
    const size_t sourceNNZ = 8;
    const std::vector<size_t> sourceRowPointers = {0, 1, 3, 5, 7, 8};
    const std::vector<size_t> sourceColIndices = {1, 0, 4, 1, 2, 0, 3, 2};
    const std::vector<double> sourceVals = {10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0};
    const CSRMatrix<double> sourceMatrix(sourceRows, sourceCols, sourceNNZ, sourceVals, sourceRowPointers, sourceColIndices);

    FixedNNZPattern thresh{2};
    SparsityPattern<double, FixedNNZPattern> pattern(sourceMatrix, thresh);
    pattern.computePattern();

    SparseApproximateMap<double, FixedNNZPattern> testSAM(targetMatrix, sourceMatrix, pattern);
    testSAM.computeMap();

    std::cout << "OK" << std::endl;
}

void testCD2D1()
{
    std::cout << "CD2D Simple Sparsity Pattern..." << std::flush;

    CSRMatrix<double> targetMatrix;
    read_mat("/home/mds222/SAM-HPC/cd2d_test/target.txt", &targetMatrix);

    CSRMatrix<double> sourceMatrix;
    read_mat("/home/mds222/SAM-HPC/cd2d_test/source.txt", &sourceMatrix);

    CSRMatrix<double> expectedMap;
    read_mat("/home/mds222/SAM-HPC/cd2d_test/cd2dtestmap.txt", &expectedMap);

    SparsityPattern<double, SimplePattern> pattern(sourceMatrix, SimplePattern());
    pattern.computePattern();
    SparseApproximateMap<double, SimplePattern> testSAM(targetMatrix, sourceMatrix, pattern);
    testSAM.computeMap();

    // assertEquals<CSRMatrix<double>>(expectedMap, *computedMap);
    std::cout << "OK" << std::endl;
}

void testCD2D2()
{
    std::cout << "CD2D Global Sparsity Pattern..." << std::flush;

    CSRMatrix<double> targetMatrix;
    read_mat("/home/mds222/SAM-HPC/cd2d_test/target.txt", &targetMatrix);

    CSRMatrix<double> sourceMatrix;
    read_mat("/home/mds222/SAM-HPC/cd2d_test/source.txt", &sourceMatrix);

    CSRMatrix<double> expectedMap;
    read_mat("/home/mds222/SAM-HPC/cd2d_test/cd2dtestmap.txt", &expectedMap);

    GlobalThresholdPattern thresh{0.001};
    SparsityPattern<double, GlobalThresholdPattern> pattern(sourceMatrix, thresh);
    pattern.computePattern();
    SparseApproximateMap<double, GlobalThresholdPattern> testSAM(targetMatrix, sourceMatrix, pattern);
    testSAM.computeMap();

    // assertEquals<CSRMatrix<double>>(expectedMap, *computedMap);
    std::cout << "OK" << std::endl;
}

void testCD2D3()
{
    std::cout << "CD2D Column Sparsity Pattern..." << std::flush;

    CSRMatrix<double> targetMatrix;
    read_mat("/home/mds222/SAM-HPC/cd2d_test/target.txt", &targetMatrix);

    CSRMatrix<double> sourceMatrix;
    read_mat("/home/mds222/SAM-HPC/cd2d_test/source.txt", &sourceMatrix);

    CSRMatrix<double> expectedMap;
    read_mat("/home/mds222/SAM-HPC/cd2d_test/cd2dtestmap.txt", &expectedMap);

    ColumnThresholdPattern thresh{0.001};
    SparsityPattern<double, ColumnThresholdPattern> pattern(sourceMatrix, thresh);
    pattern.computePattern();
    SparseApproximateMap<double, ColumnThresholdPattern> testSAM(targetMatrix, sourceMatrix, pattern);
    testSAM.computeMap();

    // assertEquals<CSRMatrix<double>>(expectedMap, *computedMap);
    std::cout << "OK" << std::endl;
}

void testCD2D4()
{
    std::cout << "CD2D Fixed NNZ Sparsity Pattern..." << std::flush;

    CSRMatrix<double> targetMatrix;
    read_mat("/home/mds222/SAM-HPC/cd2d_test/target.txt", &targetMatrix);

    CSRMatrix<double> sourceMatrix;
    read_mat("/home/mds222/SAM-HPC/cd2d_test/source.txt", &sourceMatrix);

    CSRMatrix<double> expectedMap;
    read_mat("/home/mds222/SAM-HPC/cd2d_test/cd2dtestmap.txt", &expectedMap);

    FixedNNZPattern thresh{3};
    SparsityPattern<double, FixedNNZPattern> pattern(sourceMatrix, thresh);
    pattern.computePattern();
    SparseApproximateMap<double, FixedNNZPattern> testSAM(targetMatrix, sourceMatrix, pattern);
    testSAM.computeMap();

    // assertEquals<CSRMatrix<double>>(expectedMap, *computedMap);
    std::cout << "OK" << std::endl;
}

void testTopOpt1()
{
    std::cout << "TopOpt Simple Sparsity Pattern...\n" << std::flush;

    CSRMatrix<double> targetMatrix;
    read_mat("/home/mds222/SAM-HPC/top_opt_matrices_small_csr/matrix_1.txt", &targetMatrix);
    std ::cout << "Target Matrix:\n" << targetMatrix << std::endl;

    CSRMatrix<double> sourceMatrix;
    read_mat("/home/mds222/SAM-HPC/top_opt_matrices_small_csr/matrix_2.txt", &sourceMatrix);
    std::cout << "Source Matrix:\n" << sourceMatrix << std::endl;

    SparsityPattern<double, SimplePattern> pattern(sourceMatrix, SimplePattern());
    pattern.computePattern();
    SparseApproximateMap<double, SimplePattern> testSAM(targetMatrix, sourceMatrix, pattern);
    testSAM.computeMap();

    std::cout << "OK" << std::endl;
}