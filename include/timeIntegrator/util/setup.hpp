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

template<typename VEC, typename Precise_type>
class TimeIntegrator
{
public:
    using Func = std::function<VEC(VEC,Precise_type)>;

    TimeIntegrator(){}
    ~TimeIntegrator() = default;

    virtual VEC solve(const VEC &U0, Func f, Precise_type T, int p_=0) const = 0;

private:
    mutable std::vector<VEC> result_;
};

} // namespace timeIntegrator


#endif