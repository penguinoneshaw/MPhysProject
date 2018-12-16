
#if !defined(PROJECT_FIT_H)
#define PROJECT_FIT_H
#include <vector>
#include "Minuit2/FCNBase.h"

std::vector<float> low_pass_filter(const std::vector<float> &vector, const size_t cutoff = 15);
std::vector<float> moving_average(const std::vector<float> &vector, const size_t period = 10);

class Chisquared : public ROOT::Minuit2::FCNBase {
  private:
    std::vector<double> depths;
    std::vector<double> speed_of_sound;
  public:
    Chisquared() {};
    ~Chisquared(){};

    double operator()(const std::vector<double> &par) const;

    double Up() const {
      return 1;
    };
};

#endif // PROJECT_FIT_H
