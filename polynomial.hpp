#
#if !defined(H_POLYNOMIAL)
#define H_POLYNOMIAL

#include <vector>
#include <algorithm>
#include <numeric>

namespace poly
{
  template <class T>
const T horners_method(const std::vector<T> &coefficients, const T &x_0)
{
  // Evaluates polynomials defined with coefficients vector little endian (i.e.
  // a_0 + x*a_1 + ...)

  return std::accumulate(coefficients.rbegin(), coefficients.rend(), 0.0, [x_0](T a, T b){
    return b + x_0 * a;
  });
}

template <class T>
const T horner2D(const std::vector<std::vector<T>> &coefficients, const T &x_0, const T &y_0)
{
  std::vector<T> result(coefficients.size(), 0);
  std::transform(coefficients.begin(), coefficients.end(), result.begin(),
                 [&x_0](std::vector<T> const &coeffs) -> T {
                   return horners_method<T>(coeffs, x_0);
                 });
  return horners_method(result, y_0);
}
}

#endif // H_POLYNOMIAL
