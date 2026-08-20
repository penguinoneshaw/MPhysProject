#pragma once

#include <type_traits>

namespace project::utils {
template <typename T>
concept arithmetic = std::is_arithmetic<T>::value;

}  // namespace project::utils
