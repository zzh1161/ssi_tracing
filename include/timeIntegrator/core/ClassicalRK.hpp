#ifndef TIMEINTEGRATOR_CLASSICAL_RK_HPP
#define TIMEINTEGRATOR_CLASSICAL_RK_HPP

#include "../util/setup.hpp"

namespace timeIntegrator
{

template<typename VEC, typename Precise_type>
class ClassicalRKsolver : public TimeIntegrator<VEC,Precise_type>
{
    typedef TimeIntegrator<VEC,Precise_type> Base;
    typedef Precise_type                     Real;

public:
    ClassicalRKsolver(int steps): steps_(steps){};
    ClassicalRKsolver(const ClassicalRKsolver &) = delete;
    ClassicalRKsolver& operator=(const ClassicalRKsolver &) = delete;
    ~ClassicalRKsolver() = default;

    VEC solve(const VEC &U0, Base::Func f, Real T, int p_=0) const
    {
#       define U Base::result_
        std::vector<VEC>().swap(U);
        Real k = 1.0*T/steps_;
        U.push_back(U0);
        for(int n=0; n<steps_; n++){
            auto y1 = f(U[n], n*k);
            auto y2 = f(U[n]+y1*(0.5*k), n*k+0.5*k);
            auto y3 = f(U[n]+y2*(0.5*k), n*k+0.5*k);
            auto y4 = f(U[n]+y3*k, n*k+k);
            U.push_back(U[n] + (y1+y2*2+y3*2+y4)*(1.0*k/6));
        }
        return *(U.end()-1);
#       undef U
    }

private:
    int steps_;
};

} // namespace timeIntegrator


#endif