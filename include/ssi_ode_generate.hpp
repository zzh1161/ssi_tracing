#ifndef SSI_ODE_GENERATE_HPP
#define SSI_ODE_GENERATE_HPP

#include <functional>
#include "tinynurbs/tinynurbs.h"

template<typename real>
std::function<glm::vec<4,real>(glm::vec<4,real>, real)>
ssi_ode_generate(const tinynurbs::RationalSurface<real> &p, const tinynurbs::RationalSurface<real> &q)
{
    return [&p,&q](glm::vec<4,real> u, real s){
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
        return std::move(glm::vec<4,real>(f1,f2,f3,f4));
    };
}

#endif