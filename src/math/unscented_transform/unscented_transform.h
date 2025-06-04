#pragma once

#include <cstddef>
#include <functional>
#include <memory>
#include <vector>

#include <Eigen/Core>

#include "src/math/eigen_types.hpp"
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
    UnscentedTransform(const double alpha, const double beta = 2.);
    
    std::vector<VectorNd<N>> GenerateSigmaPoints(const VectorNd<N>& state, const MatrixNd<N>& covariance);

    void Process(const VectorNd<N>& state, const MatrixNd<N>& covariance, const std::unique_ptr<Integrator>& integrator, VectorNd<N>* updated_mean, MatrixNd<N>* update_covariance);

    /**
     *
     * Getters
     *
     **/
    double alpha() const { return alpha_; }
    double beta() const { return beta_; }
    double kappa() const { return kappa_; }
    double lambda() const { return lambda_; }
    const std::vector<double>& mean_weights() const { return mean_weights_; }
    const std::vector<double>& cov_weights()  const { return cov_weights_; }

private:
    const double alpha_ = 0.;
    const double beta_ = 0.;

    double kappa_ = 0.;
    double lambda_ = 0.;

    std::vector<double> mean_weights_;
    std::vector<double> cov_weights_;
};


template<std::size_t N>
UnscentedTransform<N>::UnscentedTransform(const double alpha, const double beta) : alpha_(alpha), beta_(beta) {
    kappa_ = 3. - static_cast<double>(N);
    lambda_ = alpha_ * alpha_ * (static_cast<double>(N) + kappa_) - static_cast<double>(N);
    mean_weights_.reserve(2 * N  + 1);
    mean_weights_.emplace_back(lambda_ / (static_cast<double>(N) + lambda_));

    const double other_weights = 1 / (2 * (N + lambda_));
    for (int i = 1; i < 2 * N; ++i) {
        mean_weights_.emplace_back(other_weights);
    }

    cov_weights_ = mean_weights_;
    cov_weights_[0] += (1. - alpha_ * alpha_ + beta_);
}

template<std::size_t N>
std::vector<Eigen::Matrix<double, N, 1>> UnscentedTransform<N>::GenerateSigmaPoints(const Eigen::Matrix<double, N, 1>& state, const Eigen::Matrix<double, N, N>& covariance) {
    const Eigen::LLT<Eigen::Matrix<double, N, N>> chol(covariance);
    const Eigen::Matrix<double, N, N>  lower_cov = chol.matrixL();

    const double scaling_factor = std::sqrt(static_cast<double>(N) + lambda_);

    std::vector<Eigen::Matrix<double, N, 1>> sigma_points;
    sigma_points.reserve(2 * N + 1);
    sigma_points.emplace_back(state);

    for (int i = 1; i <= N; ++i) {
        sigma_points.emplace_back(state + scaling_factor * lower_cov.col(i));
        sigma_points.emplace_back(state - scaling_factor * lower_cov.col(i));
    }

    return sigma_points;
}

template<std::size_t N>
void Process(const VectorNd<N>& state, const MatrixNd<N>& covariance, const std::unique_ptr<Integrator>& integrator, const double timestep, VectorNd<N>* updated_mean, MatrixNd<N>* updated_covariance) {
    if (!updated_mean || !updated_covariance) {
        // TODO: Perhaps this should throw instead.
        return;
    }

    // Generate sigma points.
    const std::vector<VectorNd<N>> sigma_points = GenerateSigmaPoints(state, covariance);

    // Process sigma points with integrator.
    std::vector<VectorNd<N>> processed_sigma_points;
    processed_sigma_points.reserve(sigma_points.size());

    for (const VectorNd<N>& sigma_point : sigma_points) {
        processed_sigma_points.emplace_back(integrator->Integrate(sigma_point, timestep));
    }

    // Calculate mean.
    updated_mean = VectorNd<N>::Zeros();
    for (std::size_t i = 0; i < processed_sigma_points.size(); ++i) {
        updated_mean += mean_weights_.at(i) * processed_sigma_points.at(i);
    }

    // Now using this updated mean calculate the new covariance.
    updated_covariance = MatrixNd<N>::Zeros();
    for (std::size_t i = 0; i < processed_sigma_points.size(); ++i) {
        const VectorNd<N> difference = (processed_sigma_points.at(i) - updated_mean);
        updated_covariance += cov_weights_.at(i) * difference * difference.transpose();
    }
}

} // namespace math
