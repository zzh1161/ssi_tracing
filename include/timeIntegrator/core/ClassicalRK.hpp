#ifndef TIMEINTEGRATOR_CLASSICAL_RK_HPP
#define TIMEINTEGRATOR_CLASSICAL_RK_HPP

#include "../util/setup.hpp"

namespace timeIntegrator
{
    template<typename VEC>
    class ClassicalRKsolver : public TimeIntegrator<VEC>
    {
    public:
        using Base = TimeIntegrator<VEC>;
        using Func = std::function<VEC(VEC,Real)>;

        ClassicalRKsolver() = default;
        ClassicalRKsolver(const ClassicalRKsolver &) = delete;
        ClassicalRKsolver& operator=(const ClassicalRKsolver &) = delete;
        ~ClassicalRKsolver() = default;

        virtual VEC solve(const VEC &U0, Func f, Real T, int steps) const;
        virtual VEC solve(const VEC &U0, Func f, Real step_size) const;
    };

} // namespace timeIntegrator

#include "ClassicalRK.inl"

#endif