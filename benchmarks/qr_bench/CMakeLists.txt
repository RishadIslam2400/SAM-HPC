add_executable(qr_bench run_qr_bench.cpp)
target_link_libraries(qr_bench PUBLIC mgs_solver householder_solver sam)
target_compile_options(qr_bench
    PRIVATE

    # Debug options
    $<$<CONFIG:Debug>:-Wall>
    $<$<CONFIG:Debug>:-Werror>
    $<$<CONFIG:Debug>:-fsanitize=address>

    # Release options
    $<$<CONFIG:Release>:-O3>
    $<$<CONFIG:Release>:-DNDEBUG>
)
target_link_options(qr_bench
    PRIVATE
    $<$<CONFIG:Debug>:-fsanitize=address>
)