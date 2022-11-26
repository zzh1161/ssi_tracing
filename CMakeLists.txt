cmake_minimum_required(VERSION 3.10)
project(SSI_test)

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Libigl
include(libigl)

# Enable the target igl::glfw
igl_include(glfw)

include_directories(include)
include_directories(/usr/include/eigen3)

add_subdirectory(test)
add_subdirectory(src)