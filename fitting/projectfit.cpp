#include "projectfit.hpp"
#include <vector>
#include <cmath>
#include "fftw3.h"
#include <numeric>

std::vector<float> low_pass_filter(const std::vector<float> &vector, const size_t cutoff)
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

  for (size_t i = 0; i < in.size() / 2 + 1; i++)
  {
    float window = 1.0 ? i < cutoff : std::sqrt(std::exp(-(i - cutoff) * (i - cutoff) / 10.0));
    fft[i][0] *= window;
    fft[i][1] *= window;
  }

  fftwf_plan backward = fftwf_plan_dft_c2r_1d(in.size(), fft, out.data(), FFTW_ESTIMATE);
  fftwf_execute(backward);
  fftwf_free(fft);

  for (size_t i = 0; i < out.size(); i++)
  {
    out[i] /= out.size();
    out[i] += average;
  }

  return out;
}

std::vector<float> moving_average(const std::vector<float> &vector, const size_t period)
{
  std::vector<float> output(vector.size() - period);
  for (size_t i = 0; i < vector.size() - period; i++)
  {
    output[i] = std::accumulate(std::begin(vector) + i, std::begin(vector) + i + period, 0) / period;
  }
  return output;
}

double Chisquared::operator() (const std::vector<double> &par) const {
  return std::accumulate(par.begin(), par.end(), 0);
}