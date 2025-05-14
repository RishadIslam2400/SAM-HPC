add_executable(csr_test run_csr_tests.cpp)
target_link_libraries(csr_test PUBLIC csr_test_utils)

target_compile_options(csr_test
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
target_link_options(csr_test
    PRIVATE
    $<$<CONFIG:Debug>:-fsanitize=address>
    $<$<CONFIG:Debug>:-fsanitize=undefined>
    $<$<CONFIG:Debug>:-fsanitize=leak>
    #$<$<CONFIG:Debug>:-fsanitize=thread>
)