#
#if !defined(H_POLYNOMIAL)
#define H_POLYNOMIAL

#include <vector>
#include <algorithm>
#include <numeric>
#include "./utils.hpp"

namespace project::poly
{

  /**
   * @brief Evaluates a 1-dimensional polynomial.
   * @remark Assumes the coefficients vector provided little endian (i.e. a_0 + x*a_1 + ...)
   *
   * @tparam T the type of the coefficients of the polynomial
   * @param coefficients the coefficients of the polynomial
   * @param x_0 the point at which the polynomial should be evaluated
   * @return const T
   */
  template <project::utils::arithmetic T>
  const T horners_method(const std::vector<T> &coefficients, const T &x_0)
  {
    return std::accumulate(coefficients.rbegin(), coefficients.rend(), 0.0, [x_0](T a, T b) {
      return b + x_0 * a;
    });
  }

  /**
   * @brief Evaluates a 2-dimensional polynomial.
   * @remark Assumes the coefficients vector provided little endian (i.e. a_0 + x*a_1 + ...)
   *
   * @tparam T the type of the coefficients of the polynomial
   * @param coefficients the coefficients of the polynomial
   * @param x_0 the x-coordinate of the point at which the polynomial should be evaluated
   * @param y_0 the y-coordinate of the point at which the polynomial should be evaluated
   * @return const T
   */
  template <project::utils::arithmetic T>
  T horner2D(const std::vector<std::vector<T>> &coefficients, const T &x_0, const T &y_0)
  {
    std::vector<T> result(coefficients.size(), 0);
    std::transform(coefficients.begin(), coefficients.end(), result.begin(),
                  [&x_0](std::vector<T> const &coeffs) -> T {
                    return horners_method<T>(coeffs, x_0);
                  });
    return horners_method(result, y_0);
  }

} // namespace poly

#endif // H_POLYNOMIAL
