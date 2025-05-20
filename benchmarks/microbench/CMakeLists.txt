add_executable(microbench run_microbench.cpp)
target_link_libraries(microbench PUBLIC sam common)
target_compile_options(microbench
    PRIVATE

    # Debug options
    $<$<CONFIG:Debug>:-Wall>
    $<$<CONFIG:Debug>:-Werror>
    $<$<CONFIG:Debug>:-fsanitize=address>

    # Release options
    $<$<CONFIG:Release>:-O3>
    $<$<CONFIG:Release>:-DNDEBUG>
)
target_link_options(microbench
    PRIVATE
    $<$<CONFIG:Debug>:-fsanitize=address>
)

add_executable(threshold_benchmark run_threshold_benchmark.cpp)
target_link_libraries(threshold_benchmark PUBLIC sam common)
target_compile_options(threshold_benchmark
    PRIVATE

    # Debug options
    $<$<CONFIG:Debug>:-Wall>
    $<$<CONFIG:Debug>:-Werror>
    $<$<CONFIG:Debug>:-fsanitize=address>

    # Release options
    $<$<CONFIG:Release>:-O3>
    $<$<CONFIG:Release>:-DNDEBUG>
)
target_link_options(threshold_benchmark
    PRIVATE
    $<$<CONFIG:Debug>:-fsanitize=address>
)

option(ENABLE_SEQUENTIAL "Enable sequential execution path" OFF)
if (ENABLE_SEQUENTIAL)
    message(STATUS "Sequential execution selected")
    target_compile_definitions(microbench PRIVATE SEQUENTIAL_ENABLED)
endif()