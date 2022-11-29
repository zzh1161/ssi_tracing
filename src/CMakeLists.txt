add_executable(bezier_intersection bezier_main.cpp)
target_link_libraries(bezier_intersection PUBLIC igl::glfw)

add_executable(nurbs1_intersection nurbs1_main.cpp)
target_link_libraries(nurbs1_intersection PUBLIC igl::glfw)

add_executable(nurbs2_intersection nurbs2_main.cpp)
target_link_libraries(nurbs2_intersection PUBLIC igl::glfw)

add_executable(nurbs3_intersection nurbs3_main.cpp)
target_link_libraries(nurbs3_intersection PUBLIC igl::glfw)

add_executable(nurbs4_intersection nurbs4_main.cpp)
target_link_libraries(nurbs4_intersection PUBLIC igl::glfw)

add_executable(transversal_intersection transversal.cpp)
target_link_libraries(transversal_intersection PUBLIC igl::glfw)