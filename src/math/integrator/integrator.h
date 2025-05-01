#pragma once 

#include <functional>

#include <Eigen/Core>

#include "src/math/eigen_types.hpp"

namespace math {

template<int N>
class Integrator{
public:
    Integrator() = default;

    /**
     * Integrator abstract class
     **/
    Integrator(const std::function<VectorNd<N>(double, const VectorNd<N>&)>& func) : func_(func) {};

    virtual VectorNd<N> Integrate(const VectorNd<N>& current_state, const double timestep) = 0;

    const std::function<VectorNd<N>(double, const VectorNd<N>&)>& Get();

private:
    const std::function<VectorNd<N>(double, const VectorNd<N>&)> func_;
};

} // namespace math
