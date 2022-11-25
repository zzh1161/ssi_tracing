#include <glm/glm.hpp>
#include <tinynurbs/tinynurbs.h>
#include <timeIntegrator/timeIntegrator.h>

#include <functional>

int main()
{
    using vec4func = std::function<glm::vec4(glm::vec4, float)>;

    tinynurbs::RationalSurface3f p;
    p.degree_u = 3;
    p.degree_v = 3;
    p.knots_u = {0, 0, 0, 0, 1, 1, 1, 1};
    p.knots_v = {0, 0, 0, 0, 1, 1, 1, 1};
    p.control_points = {4, 4, 
                        {glm::vec3(-1.5,-1.5, 1.0), glm::vec3(-1.5,-0.5,1.0), glm::vec3(-1.5,0.5,1.0), glm::vec3(-1.5,1.5,1.0),
						 glm::vec3(-0.5,-1.5,1.0),  glm::vec3(-0.5,-0.5,2.0), glm::vec3(-0.5,0.5,2.0), glm::vec3(-0.5,1.5,1.0),
						 glm::vec3(0.5,-1.5,1.0),   glm::vec3(0.5,-0.5,2.0),  glm::vec3(0.5,0.5,2.0),  glm::vec3(0.5,1.5,1.0),
						 glm::vec3(1.5,-1.5,1.0),   glm::vec3(1.5,-0.5,1.0),  glm::vec3(1.5,0.5,1.0),  glm::vec3(1.5,1.5,1.0)}};
    p.weights = {4, 4,
                {1, 1, 1, 1,
                 1, 1, 1, 1,
                 1, 1, 1, 1,
                 1, 1, 1, 1}};

    tinynurbs::RationalSurface3f q;
    q.degree_u = 3;
    q.degree_v = 3;
    q.knots_u = {0, 0, 0, 0, 1, 1, 1, 1};
    q.knots_v = {0, 0, 0, 0, 1, 1, 1, 1};
    q.control_points = {4, 4, 
                        {glm::vec3(-1.5,0.5,-1.0), glm::vec3(-1.5,0.5,0.0),  glm::vec3(-1.5,0.5,1.0),  glm::vec3(-1.5,0.5,2.0),
						 glm::vec3(-0.5,0.5,-1.0), glm::vec3(-0.5,-0.5,0.0), glm::vec3(-0.5,-0.5,1.0), glm::vec3(-0.5,0.5,2.0),
					     glm::vec3(0.5,0.5,-1.0),  glm::vec3(0.5,-0.5,0.0),  glm::vec3(0.5,-0.5,1.0),  glm::vec3(0.5,0.5,2.0),
					     glm::vec3(1.5,0.5,-1.0),  glm::vec3(1.5,0.5,0.0),   glm::vec3(1.5,0.5,1.0),   glm::vec3(1.5,0.5,2.0)}};
    q.weights = {4, 4,
                {1, 1, 1, 1,
                 1, 1, 1, 1,
                 1, 1, 1, 1,
                 1, 1, 1, 1}};
    
    // sigma = u1, t = u2, u = u3, v = u4
    vec4func ode = [&p, &q](glm::vec4 u, float s){
        auto dp = tinynurbs::surfaceDerivatives(p, 1, u[0], u[1]);
        auto dpdsigma = dp(1,0), dpdt = dp(0,1);
        auto dpdu1_x_dpdu2 = glm::cross(dpdsigma, dpdt);
        auto dq = tinynurbs::surfaceDerivatives(q, 1, u[2], u[3]);
        auto dqdu = dq(1,0), dqdv = dq(0,1);
        auto dqdu3_x_dqdu4 = glm::cross(dqdu, dqdv);
        auto dcds = glm::cross(dpdu1_x_dpdu2, dqdu3_x_dqdu4);
        glm::normalize(dcds); // dc/ds
        auto module_p_sigma_t = glm::dot(dpdu1_x_dpdu2, dpdu1_x_dpdu2);
        auto module_p_u_v     = glm::dot(dqdu3_x_dqdu4, dqdu3_x_dqdu4);
        // sigma'= f1, t'= f2, u'= f3, v'= f4
        auto f1 = glm::determinant(glm::mat3(dcds, dpdt, dpdu1_x_dpdu2)) / module_p_sigma_t;
        auto f2 = glm::determinant(glm::mat3(dpdsigma, dcds, dpdu1_x_dpdu2)) / module_p_sigma_t;
        auto f3 = glm::determinant(glm::mat3(dcds, dqdv, dqdu3_x_dqdu4)) / module_p_u_v;
        auto f4 = glm::determinant(glm::mat3(dqdu, dcds, dqdu3_x_dqdu4)) / module_p_u_v;
        return std::move(glm::vec4(f1,f2,f3,f4));
    };

    timeIntegrator::ClassicalRKsolver<glm::vec4> cRKsolver;
    glm::vec4 U0{0, static_cast<float>(2.0/3.0), 0, static_cast<float>(2.0/3.0)};
    cRKsolver.solve(U0, ode, 0.0001);
    auto &U = cRKsolver.result_;
    for(int i=0; i<U.size(); ++i){
        auto P_i = tinynurbs::surfacePoint(p, U[i][0], U[i][1]);
        std::cout << P_i[0] << " " << P_i[1] << " " << P_i[2] << " ;  ";
        auto Q_i = tinynurbs::surfacePoint(q, U[i][2], U[i][3]);
        std::cout << Q_i[0] << " " << Q_i[1] << " " << Q_i[2] << std::endl;
    }

    return 0;
}