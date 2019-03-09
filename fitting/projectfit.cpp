#include "projectfit.hpp"
#include "fftw3.h"
#include "polynomial.hpp"
#include <cmath>
#include <cstddef>
#include <numeric>
#include <vector>
#include <algorithm>

#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnUserParameters.h"

namespace fit
{
std::vector<float> low_pass_filter(const std::vector<float> &vector,
                                   const std::size_t cutoff)
{
  std::vector<float> in(vector);
  float average = std::accumulate(in.begin(), in.end(), 0) / (float_t) in.size();
  for (auto &i : in)
    i -= average;
  const std::size_t FFT_ARRAY_SIZE = (in.size() / 2 + 1);
  std::vector<float> out(in.size());
  fftwf_complex *fft =
      (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * FFT_ARRAY_SIZE);
  fftwf_plan forward =
      fftwf_plan_dft_r2c_1d(in.size(), in.data(), fft, FFTW_ESTIMATE);
  fftwf_execute(forward);
  fftwf_destroy_plan(forward);

  for (std::size_t i = 0; i < FFT_ARRAY_SIZE; i++)
  {
    float window =
        1.0 ? i < cutoff
            : std::sqrt(std::exp(-(i - cutoff) * (i - cutoff) / 10.0));
    fft[i][0] *= window;
    fft[i][1] *= window;
  }

  fftwf_plan backward =
      fftwf_plan_dft_c2r_1d(in.size(), fft, out.data(), FFTW_ESTIMATE);
  fftwf_execute(backward);
  fftwf_free(fft);

  for (std::size_t i = 0; i < out.size(); i++)
  {
    out[i] /= out.size();
    out[i] += average;
  }

  return out;
}

std::vector<double_t> low_pass_filter(const std::vector<double_t> &vector,
                                      const std::size_t cutoff)
{
  std::vector<double_t> in(vector);
  double average = std::accumulate(in.begin(), in.end(), 0) / in.size();
  for (auto &i : in)
    i -= average;

  std::vector<double_t> out(in.size());
  fftw_complex *fft =
      (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * in.size() / 2 + 1);
  fftw_plan forward =
      fftw_plan_dft_r2c_1d(in.size(), in.data(), fft, FFTW_ESTIMATE);
  fftw_execute(forward);
  fftw_destroy_plan(forward);

  for (std::size_t i = 0; i < in.size() / 2 + 1; i++)
  {
    double window =
        1.0 ? i < cutoff
            : std::sqrt(std::exp(-(i - cutoff) * (i - cutoff) / 10.0));
    fft[i][0] *= window;
    fft[i][1] *= window;
  }

  fftw_plan backward =
      fftw_plan_dft_c2r_1d(in.size(), fft, out.data(), FFTW_ESTIMATE);
  fftw_execute(backward);
  fftw_free(fft);

  for (std::size_t i = 0; i < out.size(); i++)
  {
    out[i] /= out.size();
    out[i] += average;
  }

  return out;
}

template <typename T>
std::tuple<std::vector<T>, std::vector<T>> moving_average(const std::vector<T> &vector,
                                                          const std::size_t period)
{
  std::vector<T> output;
  output.reserve(vector.size());
  std::vector<T> errors;
  errors.reserve(vector.size());
  for (auto it = std::begin(vector); it + period < std::end(vector); ++it)
  {
    output.push_back(std::accumulate(it,
                                     it + period, (T)0) /
                     static_cast<T>(period));
    errors.push_back(std::accumulate(it,
                                     it + period, (T)0, [&output](auto a, auto b) { return a + (b - output.back()) * (b - output.back()); }) /
                     (T)(period - 1));
  }
  output.shrink_to_fit();
  errors.shrink_to_fit();
  return std::tie(output, errors);
}

template std::tuple<std::vector<double_t>, std::vector<double_t>>
moving_average(const std::vector<double_t> &vector, const std::size_t period);

template std::tuple<std::vector<float_t>, std::vector<float_t>> moving_average(const std::vector<float_t> &vector, const std::size_t period);

template <typename T>
T find_SOFAR_channel(const std::vector<T> &speed_of_sound, const std::vector<T> &depths, std::size_t averaging_granularity)
{
  auto [avg_sos, sos_errors] = moving_average(speed_of_sound, averaging_granularity);
  auto [avg_depths, depths_errors] = moving_average(depths, averaging_granularity);

  auto differentiate = [](const std::vector<T> &xs, const std::vector<T> &ys) {
    std::vector<T> result(ys.size(), 0), dxs(xs.size(), 0);
    std::adjacent_difference(ys.begin(), ys.end(), result.begin());
    std::adjacent_difference(xs.begin(), xs.end(), dxs.begin());
    for (size_t i = 0; i < result.size() && i < dxs.size(); i++)
    {
      result[i] = result[i] / dxs[i];
    }

    return result;
  };

  auto diff_avg_sos = differentiate(avg_depths, avg_sos);

  std::vector<size_t> maxima{0};
  if (diff_avg_sos.size() < 5)
  {
    throw std::runtime_error("NOT ENOUGH DATA");
  }

  for (auto it = diff_avg_sos.begin() + 1; it != diff_avg_sos.end(); ++it)
  {
    if (*(it - 1) > 0 && *it <= 0)
    {
      maxima.push_back(std::distance(diff_avg_sos.begin(), it));
    }
  }

  std::size_t endindex = maxima.back() == avg_sos.size() - 1 ? 0 : maxima.back();

  auto chisquared_depths = std::vector<double_t>(std::begin(avg_depths) + endindex, avg_depths.end());
  auto chisquared_sos = std::vector<double_t>(std::begin(avg_sos) + endindex, avg_sos.end());
  auto chisquared_errors = std::vector<double_t>(std::begin(sos_errors) + endindex, sos_errors.end());
  Chisquared fcn(chisquared_depths, chisquared_sos, chisquared_errors);
  ROOT::Minuit2::MnUserParameters upar;
  upar.Add("c_0", depths[0], 1000);
  upar.Add("c_1", 0, 1);
  upar.Add("c_2", 0, 1);

  ROOT::Minuit2::MnMigrad migrad(fcn, upar);
  ROOT::Minuit2::FunctionMinimum min = migrad();
  auto xmin = fcn.function_minimum(min);

  //if (!min.IsValid() || xmin > *(avg_depths.end()) || xmin < *(avg_depths.begin()) || std::isnan(xmin))

  if (!min.IsValid() || std::isnan(xmin) || xmin > avg_depths.back() || xmin < avg_depths.front())
  {
    throw std::runtime_error("out of region");
  }
  else
  {
    return (T) xmin;
  }
}

template double_t find_SOFAR_channel(const std::vector<double_t> &speed_of_sound, const std::vector<double_t> &depths, std::size_t averaging_granularity);

} // namespace fit

double bilinear_fit(const std::vector<double> &par, double d)
{
  // This is the fitting function, which I can swap out as necessary
  // Uses a bilinear fit, where the curve is theoretically fitted to two
  // successive lines, which can be generalised. Expects par to be in the format
  // (x, y, m_1, m_2) where [x, y] is the 'position' of the turning point in the
  // depth-speed plain

  double r = 0;
  if (d < par[0])
  {
    r = par[2] * (d - par[0]) + par[1];
  }
  else
  {
    r = par[3] * (d - par[0]) + par[1];
  }
  return r;
}

double polynomial_fit(const std::vector<double> &par, double d)
{
  // This is the fitting function, which I can swap out as necessary
  // Uses a cubic fitting function (originally was planned to be quadratic but
  // this way it has the same number of coefficients as the bilinear). Expects
  // par to be in the format (c_i) for i in [0, 3] is the 'position' of the
  // turning point in the depth-speed plane

  return poly::horners_method(par, d);
}

double Chisquared::operator()(const std::vector<double> &par) const
{
  std::vector<double> fitted_speeds(this->depths.size());
  std::transform(this->depths.begin(), this->depths.end(),
                 fitted_speeds.begin(),
                 [&par](double d) { return polynomial_fit(par, d); });

  double result = 0;

  for (std::size_t i = 0; i < depths.size(); i++)
  {
    result += std::pow(fitted_speeds[i] - speed_of_sound[i], 2) / this->errors[i];
  }
  return result;
}

std::vector<double>
Chisquared::fitted_to_minimisation(ROOT::Minuit2::FunctionMinimum min)
{
  std::vector<double> result(this->depths);
  auto par = min.UserState().Params();
  std::transform(result.begin(), result.end(), result.begin(),
                 [&par](double d) { return polynomial_fit(par, d); });
  return result;
}

double Chisquared::function_minimum(ROOT::Minuit2::FunctionMinimum min){
  auto minCoeff = min.UserParameters().Params();
  return -minCoeff[1] / (2 * minCoeff[2]);
}
