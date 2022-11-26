add_executable(test_glm test_glm.cpp)
add_executable(test_tn test_tn.cpp)
add_executable(test_igl test_igl.cpp)
target_link_libraries(test_igl PUBLIC igl::glfw)