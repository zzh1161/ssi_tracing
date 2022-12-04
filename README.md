<h1 align="center">Curve tracing for SSI</h1>

## File Description

+ `src` includes all SSI tracing test samples.

## Build

+ type these commands to build all samples:
    ```
    mkdir build && cd build
    cmake ..
    make
    ```

## Dependencies

+ this project uses `FetchContent` of CMake to get dependencies from [libigl](https://github.com/libigl/libigl), including `libigl`, `eigen3`, `glfw` and so on. So please **check your network** before building the project.