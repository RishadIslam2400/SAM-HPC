#include "cases/test_pattern.cpp"
#include "cases/test_sam_constructor.cpp"


int main()
{
    testSimpleSparsityPattern();
    testGlobalSparsityPattern();
    testColumnSparsityPattern();
    testFixedNNZSparsityPattern();
    testSamConstructor1();
    testSamConstructor2();
    return 0;
}