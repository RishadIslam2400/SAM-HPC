project(${CMAKE_PROJECT_NAME})
add_library(config INTERFACE config.hpp)
target_include_directories(config INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/../config/)