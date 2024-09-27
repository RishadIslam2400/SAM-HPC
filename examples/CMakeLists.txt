project(${CMAKE_PROJECT_NAME})

add_executable(sam_test SAM_test.cpp)
target_link_libraries(sam_test PUBLIC mat mx)
target_link_libraries(sam_test PUBLIC sam_func)

add_executable(householderQR_test householderQR_test.cpp)
target_link_libraries(householderQR_test PUBLIC householderQR -fsanitize=address)
set_target_properties(householderQR_test PROPERTIES COMPILE_FLAGS "-Wall -Wextra ")

add_executable(mgsQR_test mgsQR_test.cpp)
target_link_libraries(mgsQR_test PUBLIC mgsQR -fsanitize=address)
set_target_properties(mgsQR_test PROPERTIES COMPILE_FLAGS "-Wall -Wextra ")