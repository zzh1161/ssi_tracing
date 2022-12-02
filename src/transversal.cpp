#include <glm/glm.hpp>
#include <tinynurbs/tinynurbs.h>
#include <timeIntegrator/timeIntegrator.h>
#include <ssi_ode_generate.hpp>
#include <igl/opengl/glfw/Viewer.h>

#include <functional>
#include <map>
#include <vector>

int main(int argc, char *argv[])
{
    if(argc != 2)
        throw std::runtime_error("Parameter list does not meet the requirements!");
    
    tinynurbs::RationalSurface3f p, q;
    p.degree_u = 3;  q.degree_u = 3;
    p.degree_v = 3;  q.degree_v = 3;
    p.knots_u = {0, 0, 0, 0, 1, 1, 1, 1}; q.knots_u = {0, 0, 0, 0, 1, 1, 1, 1};
    p.knots_v = {0, 0, 0, 0, 1, 1, 1, 1}; q.knots_v = {0, 0, 0, 0, 1, 1, 1, 1};
    glm::vec4 U0;
    
    if(std::string(argv[1]) == "bezier"){
        p.control_points = {4, 4, 
                            {glm::vec3(-1.5,-1.5, 1.0), glm::vec3(-1.5,-0.5,1.0), glm::vec3(-1.5,0.5,1.0), glm::vec3(-1.5,1.5,1.0),
                             glm::vec3(-0.5,-1.5,1.0),  glm::vec3(-0.5,-0.5,2.0), glm::vec3(-0.5,0.5,2.0), glm::vec3(-0.5,1.5,1.0),
                             glm::vec3(0.5,-1.5,1.0),   glm::vec3(0.5,-0.5,2.0),  glm::vec3(0.5,0.5,2.0),  glm::vec3(0.5,1.5,1.0),
                             glm::vec3(1.5,-1.5,1.0),   glm::vec3(1.5,-0.5,1.0),  glm::vec3(1.5,0.5,1.0),  glm::vec3(1.5,1.5,1.0)}};
        p.weights = {4, 4, {1,1,1,1, 1,1,1,1, 1,1,1,1, 1,1,1,1}};
        q.control_points = {4, 4, 
                            {glm::vec3(-1.5,0.5,-1.0), glm::vec3(-1.5,0.5,0.0),  glm::vec3(-1.5,0.5,1.0),  glm::vec3(-1.5,0.5,2.0),
                             glm::vec3(-0.5,0.5,-1.0), glm::vec3(-0.5,-0.5,0.0), glm::vec3(-0.5,-0.5,1.0), glm::vec3(-0.5,0.5,2.0),
                             glm::vec3(0.5,0.5,-1.0),  glm::vec3(0.5,-0.5,0.0),  glm::vec3(0.5,-0.5,1.0),  glm::vec3(0.5,0.5,2.0),
                             glm::vec3(1.5,0.5,-1.0),  glm::vec3(1.5,0.5,0.0),   glm::vec3(1.5,0.5,1.0),   glm::vec3(1.5,0.5,2.0)}};
        q.weights = {4, 4, {1,1,1,1, 1,1,1,1, 1,1,1,1, 1,1,1,1}};
        U0 = {0, static_cast<float>(2.0/3.0), 0, static_cast<float>(2.0/3.0)};
    }else if(std::string(argv[1]) == "nurbs1"){
        p.control_points = {4, 4, 
                            {glm::vec3(-1.5,-1.5, 1.0), glm::vec3(-1.5,-0.5,1.0), glm::vec3(-1.5,0.5,1.0), glm::vec3(-1.5,1.5,1.0),
                             glm::vec3(-0.5,-1.5,1.0),  glm::vec3(-0.5,-0.5,2.0), glm::vec3(-0.5,0.5,2.0), glm::vec3(-0.5,1.5,1.0),
                             glm::vec3(0.5,-1.5,1.0),   glm::vec3(0.5,-0.5,2.0),  glm::vec3(0.5,0.5,2.0),  glm::vec3(0.5,1.5,1.0),
                             glm::vec3(1.5,-1.5,1.0),   glm::vec3(1.5,-0.5,1.0),  glm::vec3(1.5,0.5,1.0),  glm::vec3(1.5,1.5,1.0)}};
        p.weights = {4, 4, {1,1,1,1, 1,4,1,1, 1,1,4,1, 1,1,1,1}};
        q.control_points = {4, 4, 
                            {glm::vec3(-1.5,-1.5,0.0), glm::vec3(-1.5,-1.5,1.0), glm::vec3(-1.5,-1.5,2.0), glm::vec3(-1.5,-1.5,3.0),
                             glm::vec3(-0.5,0.5,0.0),  glm::vec3(-0.5,0.5,1.0),  glm::vec3(-0.5,0.5,2.0),  glm::vec3(-0.5,0.5,3.0),
                             glm::vec3(0.5,0.5,0.0),   glm::vec3(0.5,0.5,1.0),   glm::vec3(0.5,0.5,2.0),   glm::vec3(0.5,0.5,3.0),
                             glm::vec3(1.5,-1.5,0.0),  glm::vec3(1.5,-1.5,1.0),  glm::vec3(1.5,-1.5,2.0),  glm::vec3(1.5,-1.5,3.0)}};
        q.weights = {4, 4, {1,1,1,1, 1,20,20,1, 1,1,1,1, 1,1,1,1}};
        U0 = {0, 0, 0, static_cast<float>(1.0/3.0)};
    }else if(std::string(argv[1]) == "nurbs2"){
        p.control_points = {4, 4, 
                            {glm::vec3(-1.5,-1.5, 1.0), glm::vec3(-1.5,-0.5,1.0), glm::vec3(-1.5,0.5,1.0), glm::vec3(-1.5,1.5,1.0),
                             glm::vec3(-0.5,-1.5,1.0),  glm::vec3(-0.5,-0.5,2.0), glm::vec3(-0.5,0.5,2.0), glm::vec3(-0.5,1.5,1.0),
                             glm::vec3(0.5,-1.5,1.0),   glm::vec3(0.5,-0.5,2.0),  glm::vec3(0.5,0.5,2.0),  glm::vec3(0.5,1.5,1.0),
                             glm::vec3(1.5,-1.5,1.0),   glm::vec3(1.5,-0.5,1.0),  glm::vec3(1.5,0.5,1.0),  glm::vec3(1.5,1.5,1.0)}};
        p.weights = {4, 4, {1,1,1,1, 1,20,1,1, 1,1,20,1, 1,1,1,1}};
        q.control_points = {4, 4, 
                            {glm::vec3(-1.5,-1.5,0.0), glm::vec3(-1.5,-1.5,1.0), glm::vec3(-1.5,-1.5,2.0), glm::vec3(-1.5,-1.5,3.0),
                             glm::vec3(-1.5,0.5,0.0),  glm::vec3(-1.5,0.5,1.0),  glm::vec3(-1.5,0.5,2.0),  glm::vec3(-1.5,0.5,3.0),
                             glm::vec3(0.5,-0.5,0.0),  glm::vec3(0.5,-0.5,1.0),  glm::vec3(0.5,-0.5,2.0),  glm::vec3(0.5,-0.5,3.0),
                             glm::vec3(-1.5,-1.5,0.0), glm::vec3(-1.5,-1.5,1.0), glm::vec3(-1.5,-1.5,2.0), glm::vec3(-1.5,-1.5,3.0)}};
        q.weights = {4, 4, {1,1,1,1, 5,5,5,5, 5,5,5,5, 1,1,1,1}};
        U0 = {0, 0, 0, static_cast<float>(1.0/3.0)};
    }else if(std::string(argv[1]) == "nurbs3"){
        p.control_points = {4, 4,
                            {glm::vec3(-1.5,-1.5, 1.0), glm::vec3(-1.5,-0.5,1.0), glm::vec3(-1.5,0.5,1.0), glm::vec3(-1.5,1.5,1.0),
                             glm::vec3(-0.5,-1.5,1.0),  glm::vec3(-0.5,-0.5,1.0), glm::vec3(-0.5,0.5,1.0), glm::vec3(-0.5,1.5,1.0),
                             glm::vec3(0.5,-1.5,1.0),   glm::vec3(0.5,-0.5,1.0),  glm::vec3(0.5,0.5,1.0),  glm::vec3(0.5,1.5,1.0),
                             glm::vec3(1.5,-1.5,1.0),   glm::vec3(1.5,-0.5,1.0),  glm::vec3(1.5,0.5,1.0),  glm::vec3(1.5,1.5,1.0)}};
        p.weights = {4, 4, {1,1,1,1, 1,1,1,1, 1,1,1,1, 1,1,1,1}};
        q.control_points = {4, 4,
                            {glm::vec3(-0.5,-0.5,0.0), glm::vec3(-0.5,-0.5,1.0), glm::vec3(-0.5,-0.5,2.0), glm::vec3(-0.5,-0.5,3.0),
                             glm::vec3(-0.5,1.5,0.0),  glm::vec3(-0.5,1.5,1.0),  glm::vec3(-0.5,1.5,2.0),  glm::vec3(-0.5,1.5,3.0),
                             glm::vec3(1.5,-0.5,0.0),  glm::vec3(1.5,-0.5,1.0),  glm::vec3(1.5,-0.5,2.0),  glm::vec3(1.5,-0.5,3.0),
                             glm::vec3(-0.5,-0.5,0.0), glm::vec3(-0.5,-0.5,1.0), glm::vec3(-0.5,-0.5,2.0), glm::vec3(-0.5,-0.5,3.0)}};
        q.weights = {4, 4, {1,1,1,1, 1,3,3,15, 1,3,3,15, 1,1,1,1}};
        U0 = {static_cast<float>(1.0/3.0), static_cast<float>(1.0/3.0), 0, static_cast<float>(1.0/3.0)};
    }else if(std::string(argv[1]) == "nurbs4"){
        p.control_points = {4, 4,
                            {glm::vec3(-1.5,-1.5, 1.0), glm::vec3(-1.5,-0.5,1.0), glm::vec3(-1.5,0.5,1.0), glm::vec3(-1.5,1.5,1.0),
                             glm::vec3(-0.5,-1.5,1.0),  glm::vec3(-0.5,-0.5,2.0), glm::vec3(-0.5,0.5,2.0), glm::vec3(-0.5,1.5,1.0),
                             glm::vec3(0.5,-1.5,1.0),   glm::vec3(0.5,-0.5,2.0),  glm::vec3(0.5,0.5,2.0),  glm::vec3(0.5,1.5,1.0),
                             glm::vec3(1.5,-1.5,1.0),   glm::vec3(1.5,-0.5,1.0),  glm::vec3(1.5,0.5,1.0),  glm::vec3(1.5,1.5,1.0)}};
        p.weights = {4, 4, {1,1,1,1, 1,1,20,1, 1,20,1,1, 1,1,1,1}};
        q.control_points = {4, 4,
                            {glm::vec3(-0.5,-2.5,0.0), glm::vec3(-0.5,-2.5,1.0), glm::vec3(-0.5,-2.5,2.0), glm::vec3(-0.5,-2.5,3.0),
                             glm::vec3(0.5,-1.5,0.0),  glm::vec3(0.5,-1.5,1.0),  glm::vec3(0.5,-1.5,2.0),  glm::vec3(0.5,-1.5,3.0),
                             glm::vec3(1.5,-0.5,0.0),  glm::vec3(1.5,-0.5,1.0),  glm::vec3(1.5,-0.5,2.0),  glm::vec3(1.5,-0.5,3.0),
                             glm::vec3(2.5,0.5,0.0),   glm::vec3(2.5,0.5,1.0),   glm::vec3(2.5,0.5,2.0),   glm::vec3(2.5,0.5,3.0)}};
        q.weights = {4, 4, {1,1,1,1, 1,1,1,1, 1,1,1,1, 1,1,1,1}};
        U0 = {static_cast<float>(2.0/3.0), 0, static_cast<float>(1.0/3.0), static_cast<float>(1.0/3.0)};
    }else{
        throw std::runtime_error("No such test sample!");
    }

    auto ode = ssi_ode_generate(p, q);
    timeIntegrator::ClassicalRKsolver<glm::vec4> cRKsolver;
    cRKsolver.solve(U0, ode, 0.0001);
    auto &U = cRKsolver.result_;
    Eigen::MatrixXd addPoints(U.size(),3);
    for(int i=0; i<U.size(); ++i){
        auto P_i = tinynurbs::surfacePoint(p, U[i][0], U[i][1]);
        // std::cout << P_i[0] << " " << P_i[1] << " " << P_i[2] << " ;  ";
        // auto Q_i = tinynurbs::surfacePoint(q, U[i][2], U[i][3]);
        // std::cout << Q_i[0] << " " << Q_i[1] << " " << Q_i[2] << std::endl;
        addPoints(i,0) = P_i[0];
        addPoints(i,1) = P_i[1];
        addPoints(i,2) = P_i[2];
    }

    igl::opengl::glfw::Viewer viewer;
    std::vector<std::string> names;
    if(std::string(argv[1]) == "bezier")
        names = {"surface1.obj", "surface2.obj"};
    else
        names = {"surface1.stl", "surface2.stl"};
    std::map<int, Eigen::RowVector3d> colors;
    int last_selected = -1;
    for(const auto & name : names){
        viewer.load_mesh_from_file(std::string("../src/examples/") + 
                                   std::string(argv[1]) + std::string("/") + name);
        colors.emplace(viewer.data().id, 0.5*Eigen::RowVector3d::Random().array()+0.5);
    }
    viewer.callback_key_down =
        [&](igl::opengl::glfw::Viewer &, unsigned int key, int mod)
    {
        if(key == GLFW_KEY_BACKSPACE){
            int old_id = viewer.data().id;
            if (viewer.erase_mesh(viewer.selected_data_index)){
                colors.erase(old_id);
                last_selected = -1;
            }
            return true;
        }
        return false;
    };
    viewer.callback_pre_draw =
        [&](igl::opengl::glfw::Viewer &)
    {
        if (last_selected != viewer.selected_data_index){
            for (auto &data : viewer.data_list){
                data.set_colors(colors[data.id]);
            }
            viewer.data_list[viewer.selected_data_index].set_colors(Eigen::RowVector3d(0.9,0.1,0.1));
            last_selected = viewer.selected_data_index;
        }
        return false;
    };
    viewer.data().add_points(addPoints, Eigen::RowVector3d(0, 1, 0));
    viewer.data().point_size = 3.0;
    viewer.launch();

    return 0;
}