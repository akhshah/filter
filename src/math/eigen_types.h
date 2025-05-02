#pragma once

#include <cstdef>
#include <vector>

#include <Eigen/Core>

namespace math {

using Vector2d = Eigen::Vector2d;
using Vector3d = Eigen::Vector3d;

template<std::size_t N>
using VectorNd = Eigen::Matrix<double, N, 1>;

template<std::size_t N>
using MatrixNd = Eigen::Matrix<double, N, N>;

} // namespace math
