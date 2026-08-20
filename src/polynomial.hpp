#pragma once

#include <algorithm>
#include <concepts>
#include <functional>
#include <numeric>
#include <ranges>
#include <type_traits>
#include <vector>

#include "./utils.hpp"
namespace project::poly {

/**
 * @brief Evaluates a 1-dimensional polynomial.
 * @remark Assumes the coefficients vector provided little endian (i.e. a_0 +
 * x*a_1 + ...)
 *
 * @tparam T the type of the coefficients of the polynomial
 * @param coefficients the coefficients of the polynomial
 * @param x_0 the point at which the polynomial should be evaluated
 * @return const T
 */
[[nodiscard]] auto horners_method(
    const std::ranges::range auto& coefficients,
    const std::ranges::range_value_t<decltype(coefficients)>& x_0) {
    using T = std::ranges::range_value_t<decltype(coefficients)>;
    T value{};

    for (auto coefficient : coefficients) {
        value += x_0 * coefficient;
    }

    return std::ranges::fold_left(coefficients | std::views::transform([x_0](auto v) { return v * x_0; }), 0, std::plus<>{});
}

/**
 * @brief Evaluates a 2-dimensional polynomial.
 * @remark Assumes the coefficients vector provided little endian (i.e. a_0 +
 * x*a_1 + ...)
 *
 * @tparam T the type of the coefficients of the polynomial
 * @param coefficients the coefficients of the polynomial
 * @param x_0 the x-coordinate of the point at which the polynomial should be
 * evaluated
 * @param y_0 the y-coordinate of the point at which the polynomial should be
 * evaluated
 * @return const T
 */
template <project::utils::arithmetic T>
T horner2D(const std::vector<std::vector<T>>& coefficients, const T& x_0,
           const T& y_0) {
    return horners_method(coefficients | std::views::transform([&x_0](std::vector<T> const& coeffs) -> T {
                              return horners_method(coeffs, x_0);
                          }),
                          y_0);
}

}  // namespace project::poly
