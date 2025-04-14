#include "cases/test_pattern.cpp"


int main()
{
    testSimpleSparsityPattern();
    testGlobalSparsityPattern();
    testColumnSparsityPattern();
    testFixedNNZSparsityPattern();
    return 0;
}