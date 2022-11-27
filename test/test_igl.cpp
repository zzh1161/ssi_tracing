#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h>

int main(int argc, char *argv[])
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    igl::readOFF("../test/examples/bunny.off", V,F);
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(V,F);
    viewer.launch();

    return 0;
}
