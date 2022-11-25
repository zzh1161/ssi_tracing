#ifndef TIMEINTEGRATOR_SETUP_HPP
#define TIMEINTEGRATOR_SETUP_HPP

#include<cassert>
#include<cmath>
#include<deque>
#include<iostream>
#include<vector>
#include<functional>
#include<memory>

namespace timeIntegrator
{
    using Real = float;

    template<typename VEC>
    class TimeIntegrator
    {
    public:
        using Func = std::function<VEC(VEC,Real)>;

        TimeIntegrator(){}
        ~TimeIntegrator() = default;

        virtual VEC solve(const VEC &U0, Func f, Real T, int steps) const = 0;
        virtual VEC solve(const VEC &U0, Func f, Real step_size) const = 0;

    public:
        mutable std::vector<VEC> result_;
    };

    /**
     * two ways to judge whether to exit tracing
     * 1. if U equals to U0, it's a loop.
     * 2. when U beyond the bound of the parameter field.
     */
    template<typename VEC>
    class ExitCondition
    {
    public:
        static constexpr Real epsilon = 0.00000001;

    public:
        ExitCondition() = default;
        ExitCondition(VEC U0): _U0(U0){}
        ~ExitCondition() = default;

        bool operator()(const VEC &U) const {
            VEC err = U-_U0;
            Real error = static_cast<Real>(std::abs(err[0]) + std::abs(err[1]) +
                                           std::abs(err[2]) + std::abs(err[3]));
            if(error < epsilon)
                return true;
            if((U[0]<=0 || U[0]>=1) || (U[1]<=0 || U[1]>=1) ||
            (U[2]<=0 || U[2]>=1) || (U[3]<=0 || U[3]>=1))
                return true;
            return false;
        }

    private:
        VEC _U0;
    };

} // namespace timeIntegrator


#endif