#pragma once

#include <functional>

#include <Eigen/Core>

#include "integrator.h"

namespace math {
namespace integrator {

template<int N>
class ForwardEuler : public Integrator<N> {
public:
    ForwardEuler() = default;
    ForwardEuler(const std::function<Eigen::Matrix<double, N, 1>(double, const Eigen::Matrix<double, N, 1>&)>& func) : Integrator(func) {};

    Eigen::Matrix<double, N, 1> Integrate(const Eigen::Matrix<double, N, 1>& curr_state, const double timestep) override {
        return curr_state + Get()(timestep, curr_state);
    }
};

} // namespace integrator
} // namespace math
