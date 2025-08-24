#include "cases/test_pattern.cpp"
#include "cases/test_read_mat.cpp"
#include "cases/test_sam_computation.cpp"
#include "cases/test_sam_solver.cpp"

int main()
{
    testSimpleSparsityPattern();
    testGlobalSparsityPattern1();
    testGlobalSparsityPattern2();
    testColumnSparsityPattern();
    testFixedNNZSparsityPattern();
    testReadMat1();
    testReadMat2();
    testSAMSanityCheck1();
    testSAMSanityCheck2();
    testSAMSanityCheck3();
    testCD2D1();
    testCD2D2();
    testCD2D3();
    testCD2D4();
    // testTopOpt1();
    test_sam_solver_cd2d_1();
    test_sam_solver_cd2d_2();
    test_sam_solver_cd2d_3();
    test_sam_solver_cd2d_4();
    test_sam_solver_top_opt_1();
    return 0;
}