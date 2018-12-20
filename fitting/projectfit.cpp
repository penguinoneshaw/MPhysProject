#include "projectfit.hpp"
#include <cstddef>
#include <vector>
#include <cmath>
#include "fftw3.h"
#include <numeric>
#include "polynomial.hpp"

std::vector<float> low_pass_filter(const std::vector<float> &vector, const std::size_t cutoff)
{
  std::vector<float> in(vector);
  float average = std::accumulate(in.begin(), in.end(), 0) / in.size();
  for (auto &i : in)
    i -= average;

  std::vector<float> out(in.size());
  fftwf_complex *fft = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * in.size() / 2 + 1);
  fftwf_plan forward = fftwf_plan_dft_r2c_1d(in.size(), in.data(), fft, FFTW_ESTIMATE);
  fftwf_execute(forward);
  fftwf_destroy_plan(forward);

  for (std::size_t i = 0; i < in.size() / 2 + 1; i++)
  {
    float window = 1.0 ? i < cutoff : std::sqrt(std::exp(-(i - cutoff) * (i - cutoff) / 10.0));
    fft[i][0] *= window;
    fft[i][1] *= window;
  }

  fftwf_plan backward = fftwf_plan_dft_c2r_1d(in.size(), fft, out.data(), FFTW_ESTIMATE);
  fftwf_execute(backward);
  fftwf_free(fft);

  for (std::size_t i = 0; i < out.size(); i++)
  {
    out[i] /= out.size();
    out[i] += average;
  }

  return out;
}

std::vector<float> moving_average(const std::vector<float> &vector, const std::size_t period)
{
  std::vector<float> output(vector.size() - period);
  for (std::size_t i = 0; i < vector.size() - period; i++)
  {
    output[i] = std::accumulate(std::begin(vector) + i, std::begin(vector) + i + period, 0) / period;
  }
  return output;
}

constexpr double bilinear_fit(const std::vector<double> &par, double d) { 
  // This is the fitting function, which I can swap out as necessary
  // Uses a bilinear fit, where the curve is theoretically fitted to two successive lines, which can be generalised.
  // Expects par to be in the format (x, y, m_1, m_2) where [x, y] is the 'position' of the turning point in the depth-speed plain
   
  double r = 0;
  if (d < par[0]){
    r = par[2] * (d - par[0]) + par[1];
  } else {
    r = par[3] * (d - par[0]) + par[1];
  }
  return r;
}

constexpr double polynomial_fit(const std::vector<double> &par, double d) { 
  // This is the fitting function, which I can swap out as necessary
  // Uses a cubic fitting function (originally was planned to be quadratic but this way it has the same number of coefficients as the bilinear).
  // Expects par to be in the format (c_i) for i in [0, 3] is the 'position' of the turning point in the depth-speed plain
   
  return poly::horners_method(par, d);
}

double Chisquared::operator() (const std::vector<double> &par) const {
  std::vector<double> fitted_speeds(this->depths.size());
  std::transform(this->depths.begin(), this->depths.end(), fitted_speeds.begin(), [&par](double d) {
    return polynomial_fit(par, d);
  });

  double result = 0;
  
  for (std::size_t i = 0; i < depths.size(); i++){
    result += std::pow(fitted_speeds[i] - speed_of_sound[i], 2);
  }
  return result;
}

std::vector<double> Chisquared::fitted_to_minimisation(ROOT::Minuit2::FunctionMinimum min){
  std::vector<double> result(this->depths);
  auto par = min.UserState().Params();
  std::transform(result.begin(), result.end(), result.begin(), [&par](double d) {
    return polynomial_fit(par, d);
  });
  return result;
};
