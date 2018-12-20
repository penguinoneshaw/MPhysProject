#
#if !defined(H_SPEED_OF_SOUND)
#define H_SPEED_OF_SOUND

#include <algorithm>
#include <cmath>
#include <vector>

#include "polynomial.hpp"

namespace speed_of_sound {
template <class T>
const std::vector<std::vector<T>> C_coeffs{
    {1402.388, 5.03830, -5.81090e-2, 3.3432e-4, -1.47797e-6, 3.1419e-9},
    {0.153563, 6.8999e-4, -8.1829e-6, 1.3632e-7, -6.1260e-10, 0},
    {3.1260e-5, -1.7111e-6, 2.5986e-8, -2.5353e-10, 1.0415e-12, 0},
    {-9.7729e-9, 3.8513e-10, -2.3654e-12, 0, 0, 0}};

template <class T>
const std::vector<std::vector<T>> A_coeffs{
    {1.389, -1.262E-2, 7.166E-5, 2.008E-6, -3.21E-8},
    {9.4742E-5, -1.2583E-5, -6.4928E-8, 1.0515E-8, -2.0142E-10},
    {-3.9064E-7, 9.1061E-9, -1.6009E-10, 7.994E-12, 0},
    {1.100E-10, 6.651E-12, -3.391E-13, 0, 0}};

template <class T>
const std::vector<std::vector<T>> B_coeffs{{-1.922E-2, -4.42E-5},
                                           {7.3637E-5, 1.7950E-7}};
                                           
template <class T> 
const T pressure_at_depth(T depth, T latitude = -30.0) {
  latitude = latitude * M_PI / 180.0;
  T h = depth * (1.00818e-2 +
                 depth * (2.465e-8 + depth * (-1.25e-13 + depth * 2.8e-19)));
  T k = (9.7803 * (1.0 + 5.3e-3 * std::pow(std::sin(latitude), 2)) -
         2e-5 * depth) /
        (9.80612 - 2e-5 * depth);
  T thyh_0 = (1.0e-2 / (depth + 100) + 6.2e-6) * depth;
  return 1000.0 * (h * k - thyh_0);
}

template <class T>
const T speed_of_sound(const T &p, const T &t, const T &s) {
  auto p_kpa = p / 100;
  using std::vector;

  const auto D = 1.727E-3 - 7.9836E-6 * p_kpa;
  const auto C_w = poly::horner2D(C_coeffs<T>, t, p_kpa);

  const auto B = poly::horner2D(B_coeffs<T>, t, p_kpa);
  const auto A = poly::horner2D(A_coeffs<T>, t, p_kpa);

  return C_w + A * s + B * s * std::sqrt(s) + D * s * s;
}
} // namespace speed_of_sound

#endif // H_SPEED_OF_SOUND
