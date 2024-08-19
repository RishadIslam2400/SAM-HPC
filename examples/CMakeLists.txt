project(${CMAKE_PROJECT_NAME})

add_executable(sam_test SAM_test.cpp)
target_link_libraries(sam_test PUBLIC mat mx)
target_link_libraries(sam_test PUBLIC sam_func)