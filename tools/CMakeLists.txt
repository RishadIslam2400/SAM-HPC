add_library(tools INTERFACE)
target_include_directories(tools 
    INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../tools/include/>
)