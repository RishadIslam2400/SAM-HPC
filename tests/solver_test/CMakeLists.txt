add_executable(solver_test run_solver_tests.cpp)
target_link_libraries(solver_test PUBLIC solver preconditioner common test_utils)
target_compile_options(solver_test
    PRIVATE

    # Debug options
    $<$<CONFIG:Debug>:-Wall>
    $<$<CONFIG:Debug>:-Werror>
    $<$<CONFIG:Debug>:-fsanitize=address>
    $<$<CONFIG:Debug>:-fsanitize=undefined>
    $<$<CONFIG:Debug>:-fsanitize=leak>
    #$<$<CONFIG:Debug>:-fsanitize=thread>
    $<$<CONFIG:Debug>:-g>
    $<$<CONFIG:Debug>:-DDEBUG>

    # Release options
    $<$<CONFIG:Release>:-O3>
    $<$<CONFIG:Release>:-DNDEBUG>
)
target_link_options(solver_test
    PRIVATE
    $<$<CONFIG:Debug>:-fsanitize=address>
    $<$<CONFIG:Debug>:-fsanitize=undefined>
    $<$<CONFIG:Debug>:-fsanitize=leak>
    #$<$<CONFIG:Debug>:-fsanitize=thread>
)