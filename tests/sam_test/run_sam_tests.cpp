#include "cases/test_pattern.cpp"
#include "cases/test_sam_constructor.cpp"
#include "cases/test_read_mat.cpp"
#include "cases/test_sam_computation.cpp"


int main()
{
    testSimpleSparsityPattern();
    testGlobalSparsityPattern1();
    testGlobalSparsityPattern2();
    testColumnSparsityPattern();
    testFixedNNZSparsityPattern();
    testSamConstructor1();
    testSamConstructor2();
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
    return 0;
}