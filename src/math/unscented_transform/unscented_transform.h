#pragma once

#include <functional>

#include <Eigen/Core>

#include "src/math/integrator/integrator.h"

namespace math {

/**
 *
 *
 **/
template<std::size_t N>
class UnscentedTransform {
public:
    UnscentedTransform() = default;
    UnscentedTransform(const double alpha, const double beta);
    
    std::vector<Eigen::Matrix<double, N, 1>> GenerateSigmaPoints(const Eigen::Matrix<double, N, 1> state, const Eigen::Matrix<double, N, N>);

    std::vector<Eigen::Matrix<double, N, 1>> Process();

private:
    const double alpha = 0.;
    const double beta = 0.;

    double kappa = 0.;
    double lambda = 0.;
};

} // namespace math

