add_executable(bezier_intersection bezier_main.cpp)
target_link_libraries(bezier_intersection PUBLIC igl::glfw)

add_executable(nurbs1_intersection nurbs1_main.cpp)
target_link_libraries(nurbs1_intersection PUBLIC igl::glfw)

add_executable(nurbs2_intersection nurbs2_main.cpp)
target_link_libraries(nurbs2_intersection PUBLIC igl::glfw)