// ClassicalRK_inl

namespace timeIntegrator
{
    template<typename VEC>
    VEC ClassicalRKsolver<VEC>::solve(const VEC &U0,
            ClassicalRKsolver<VEC>::Func f, Real T, int steps) const
    {
        auto &U = ClassicalRKsolver<VEC>::Base::result_;
        U.clear();
        Real k = 1.0*T/steps;
        U.push_back(U0);
        for(int n=0; n<steps; n++){
            auto y1 = f(U[n], n*k);
            auto y2 = f(U[n]+y1*static_cast<Real>(0.5*k), n*k+0.5*k);
            auto y3 = f(U[n]+y2*static_cast<Real>(0.5*k), n*k+0.5*k);
            auto y4 = f(U[n]+y3*k, n*k+k);
            U.push_back(U[n] + (y1+y2*static_cast<Real>(2)+y3*static_cast<Real>(2)+y4)*static_cast<Real>(1.0*k/6));
        }
        return U.back();
    }

    template<typename VEC>
    VEC ClassicalRKsolver<VEC>::solve(const VEC &U0, 
            ClassicalRKsolver<VEC>::Func f, Real step_size) const
    {
        auto &U = ClassicalRKsolver<VEC>::Base::result_;
        U.clear();
        ExitCondition<VEC> exit_c(U0);
        Real k      = step_size;
        auto next_U = U0;
        int  n      = 0;
        do{
            U.push_back(next_U);
            auto y1 = f(U[n], n*k);
            auto y2 = f(U[n]+y1*static_cast<Real>(0.5*k), n*k+0.5*k);
            auto y3 = f(U[n]+y2*static_cast<Real>(0.5*k), n*k+0.5*k);
            auto y4 = f(U[n]+y3*k, n*k+k);
            next_U = U[n] + (y1+y2*static_cast<Real>(2)+y3*static_cast<Real>(2)+y4)*static_cast<Real>(1.0*k/6);
            n++;
        }while(!exit_c(next_U));
        return U.back();
    }

} // namespace timeIntegrator
