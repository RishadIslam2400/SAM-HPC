project(${CMAKE_PROJECT_NAME})

add_executable(sam_test SAM_test.cpp)
target_link_libraries(sam_test PUBLIC mat mx)
target_link_libraries(sam_test PUBLIC sam_func)

add_executable(sam_test_asan SAM_test.cpp)
target_link_libraries(sam_test_asan PUBLIC mat mx -fsanitize=address)
target_link_libraries(sam_test_asan PUBLIC sam_func)

add_executable(householderQR_test householderQR_test.cpp)
target_link_libraries(householderQR_test PUBLIC householderQR -fsanitize=address)
set_target_properties(householderQR_test PROPERTIES COMPILE_FLAGS "-Wall -Wextra ")

add_executable(mgsQR_test mgsQR_test.cpp)
target_link_libraries(mgsQR_test PUBLIC mgsQR -fsanitize=address)
set_target_properties(mgsQR_test PROPERTIES COMPILE_FLAGS "-Wall -Wextra ")

add_executable(testSparstiyPattern testSparsityPattern.cpp)
target_link_libraries(testSparstiyPattern PUBLIC mat mx sam_func)
set_target_properties(testSparstiyPattern PROPERTIES COMPILE_FLAGS "-Wall -Wextra ")

add_executable(test_spgemm test_spgemm.cpp)
target_link_libraries(test_spgemm PUBLIC mat mx sam_func -fsanitize=address)
set_target_properties(test_spgemm PROPERTIES COMPILE_FLAGS "-Wall -Wextra ")