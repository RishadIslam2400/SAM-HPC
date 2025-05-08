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

    assertEquals<std::string>("(0, 0): 1\n(0, 2): 4\n(0, 3): 5\n(1, 0): 2\n(1, 1): -1\n(2, 2): 3\n(2, 3): 2\n", oss.str());

    m.set(7, 0, 2).set(5, 1, 1).set(3, 2, 0);

    oss.str("");
    oss << m;
    assertEquals<std::string>("(0, 0): 1\n(0, 2): 7\n(0, 3): 5\n(1, 0): 2\n(1, 1): 5\n(2, 0): 3\n(2, 2): 3\n(2, 3): 2\n", oss.str());

    std::cout << " OK" << std::endl;
}