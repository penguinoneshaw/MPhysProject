#if !defined(H_PROJECT_UTILS)
#define H_PROJECT_UTILS

#include <type_traits>

namespace project::utils
{
    template <typename T>
    concept arithmetic = std::is_arithmetic<T>::value;

} // namespace project::utils

#endif // !defined(H_PROJECT_UTILS)