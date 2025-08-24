add_subdirectory(qr_bench)
add_subdirectory(microbench)

add_executable(top_opt_bench topology_benchmark.cpp)
target_link_libraries(top_opt_bench PUBLIC sam common solver)
target_compile_options(top_opt_bench
    PRIVATE

    # Debug options
    $<$<CONFIG:Debug>:-Wall>
    $<$<CONFIG:Debug>:-Werror>
    $<$<CONFIG:Debug>:-fsanitize=address>

    # Release options
    $<$<CONFIG:Release>:-O3>
    $<$<CONFIG:Release>:-DNDEBUG>
)
target_link_options(top_opt_bench
    PRIVATE
    $<$<CONFIG:Debug>:-fsanitize=address>
)

add_executable(weather_bench weather_benchmark.cpp)
target_link_libraries(weather_bench PUBLIC sam common solver)
target_compile_options(weather_bench
    PRIVATE

    # Debug options
    $<$<CONFIG:Debug>:-Wall>
    $<$<CONFIG:Debug>:-Werror>
    $<$<CONFIG:Debug>:-fsanitize=address>

    # Release options
    $<$<CONFIG:Release>:-O3>
    $<$<CONFIG:Release>:-DNDEBUG>
)
target_link_options(weather_bench
    PRIVATE
    $<$<CONFIG:Debug>:-fsanitize=address>
)