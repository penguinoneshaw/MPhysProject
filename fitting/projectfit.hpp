
#if !defined(PROJECT_FIT_H)
#define PROJECT_FIT_H
#include <numeric>
#include <vector>
#include <cstddef>
#include "Minuit2/FCNBase.h"
#include "Minuit2/FunctionMinimum.h"
#include "fftw3.h"
#include <map>
#include <tuple>

namespace fit
{
std::vector<float> low_pass_filter(const std::vector<float> &vector, const std::size_t cutoff = 15);
template <typename T>
std::tuple<std::vector<T>, std::vector<T>> moving_average(const std::vector<T> &vector, const std::size_t period = 10);
template <typename T>
T find_SOFAR_channel(const std::vector<T> &speed_of_sound, const std::vector<T> &depths);

template <typename K, typename V>
std::vector<V> analyse_periodicity(std::map<K, V> t_series_data)
{
  // Assumes integer spaced time series
  std::vector<V> input_vector{};
  input_vector.reserve(25535);
  V avg = std::accumulate(t_series_data.begin(), t_series_data.end(), 0, [](auto a, auto b){return a + b.second;}) / (V) t_series_data.size();
  V var = std::accumulate(t_series_data.begin(), t_series_data.end(), 0, [avg](auto a, auto b){return a + (b.second - avg)*(b.second - avg);})/(V) (t_series_data.size() - 1);
  auto std_dev = std::sqrt(var);
  auto it = t_series_data.begin();
  std::pair<K, V> prev = *it;
  input_vector.push_back((prev.second - avg)/std_dev);
  for (std::next(it); it != std::end(t_series_data); std::next(it)) {
    std::size_t abs_diff = std::abs((int64_t) it->first - (int64_t) prev.first);
    if (abs_diff > 1) {
      for (std::size_t i = 0; i < abs_diff; i++) {
        input_vector.push_back(0);
      }
    }
    prev = std::pair<K, V>(it->first, it->second);
    input_vector.push_back((prev.second - avg)/std_dev);
  }

  return input_vector;
}
} // namespace fit

class Chisquared : public ROOT::Minuit2::FCNBase
{
private:
  std::vector<double> depths;
  std::vector<double> speed_of_sound;
  std::vector<double> errors;

public:
  Chisquared(std::vector<double> depths, std::vector<double> speed_of_sound, std::vector<double> errors) : depths{depths}, speed_of_sound{speed_of_sound}, errors{errors} {};
  ~Chisquared(){};

  double operator()(const std::vector<double> &par) const;

  std::vector<double> fitted_to_minimisation(ROOT::Minuit2::FunctionMinimum min);

  double Up() const
  {
    return 1;
  };
};

#endif // PROJECT_FIT_H
