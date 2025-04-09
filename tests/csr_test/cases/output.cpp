#include "CSRMatrixMock.hpp"
#include "testlib.hpp"

void testOutput()
{
    std::cout << "Output..." << std::flush;

    std::ostringstream oss;
    std::string output;

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
    CSRMatrixMock<int> m(3, 4, vals1, rowPointers1, colIndices1);
    
    oss << m;

    assertEquals<std::string>("1 0 4 5\n2 -1 0 0\n0 0 3 2", oss.str());

    m.set(7, 0, 2).set(5, 1, 1).set(3, 2, 0);

    oss.str("");
    oss << m;
    assertEquals<std::string>("1 0 7 5\n2 5 0 0\n3 0 3 2", oss.str());

    std::cout << " OK" << std::endl;
}