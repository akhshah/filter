#pragma once 

#include <functional>

#include <Eigen/Core>

namespace math {
namespace integrator {

template<int N>
class Integrator{
public:
    Integrator() = default;

    /**
     * Integrator abstract class
     **/
    Integrator(const std::function<Eigen::Matrix<double, N, 1>(double, const Eigen::Matrix<double, N, 1>&)>& func) : func_(func) {};

    virtual Eigen::Matrix<double, N, 1> Integrate(const Eigen::Matrix<double, N, 1>& current_state, const double timestep) = 0;

    const std::function<Eigen::Matrix<double, N, 1>(double, const Eigen::Matrix<double, N, 1>&)>& Get();

private:
    const std::function<Eigen::Matrix<double, N, 1>(double, const Eigen::Matrix<double, N, 1>&)> func_;
};

} // namespace integrator
} // namespace math
