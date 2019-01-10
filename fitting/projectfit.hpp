
#if !defined(PROJECT_FIT_H)
#define PROJECT_FIT_H
#include <vector>
#include <cstddef>
#include "Minuit2/FCNBase.h"
#include "Minuit2/FunctionMinimum.h"

namespace fit {
  std::vector<float> low_pass_filter(const std::vector<float> &vector, const std::size_t cutoff = 15);
  template <typename T> std::vector<T> moving_average(const std::vector<T> &vector, const std::size_t period = 10);
  template <typename T> std::size_t find_SOFAR_channel(const std::vector<T> &speed_of_sound);
}

class Chisquared : public ROOT::Minuit2::FCNBase
{
private:
  std::vector<double> depths;
  std::vector<double> speed_of_sound;

public:
  Chisquared(std::vector<double> depths, std::vector<double> speed_of_sound) : depths{depths}, speed_of_sound{speed_of_sound} {};
  ~Chisquared(){};

  double operator()(const std::vector<double> &par) const;

  std::vector<double> fitted_to_minimisation(ROOT::Minuit2::FunctionMinimum min);

  double Up() const
  {
    return 1;
  };
};

#endif // PROJECT_FIT_H
