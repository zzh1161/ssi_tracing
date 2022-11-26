# Curve tracing for SSI

## File Description

+ `src` includes the first SSI tracing test.

## Build

+ type these commands to build all examples:
    ```
    mkdir build && cmake ..
    make
    ```

## Dependencies

+ this project uses `FetchContent` of CMake to get dependencies from [libigl](https://github.com/libigl/libigl), including `libigl`, `eigen3`, `glfw` and so on. So please **check your network** before building the project.