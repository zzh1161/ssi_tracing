cmake_minimum_required(VERSION 3.10)
project(SSI_test)

set(CMAKE_CXX_STANDARD 17)

include_directories(include)
include_directories(/usr/include/eigen3)
include_directories(/usr/loacl/include)

add_subdirectory(test)
add_subdirectory(src)