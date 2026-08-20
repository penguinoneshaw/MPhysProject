#pragma once

#include <cmath>
#include <complex>
#include <cstddef>
#include <map>
#include <numeric>
#include <tuple>
#include <utility>
#include <vector>

#include "Minuit2/FCNBase.h"
#include "Minuit2/FunctionMinimum.h"
#include "fftw3.h"

namespace project {

namespace fit {

enum FitFunction { FIT_QUADRATIC,
                   FIT_IDEAL };

std::vector<float> low_pass_filter(const std::vector<float>& vector,
                                   const std::size_t cutoff = 15);
std::vector<double> low_pass_filter(const std::vector<double>& vector,
                                    const std::size_t cutoff = 15);

double ideal_sound_channel(const std::vector<double>& par, double z);

template <FitFunction f>
ROOT::Minuit2::MnUserParameters
build_user_parameters(ROOT::Minuit2::MnUserParameters& user_parameters,
                      double estimated_minimum_depth = 500,
                      double estimated_minimum_speed = 1500) {
    throw std::runtime_error("No specified algorithm");
}

template <typename T, FitFunction f = FIT_QUADRATIC>
std::pair<T, T> find_sofar_channel(const std::vector<T>& speed_of_sound,
                                   const std::vector<T>& depths,
                                   std::size_t averaging_granularity = 10);

template <typename K, typename V>
std::map<double, std::complex<double>>
analyse_periodicity(std::map<K, V> t_series_data, std::size_t block_size = 1) {
    // Assumes integer spaced time series
    std::vector<V> input_vector{};
    if (t_series_data.size() < 10) {
        throw std::runtime_error("Not enough data to FFT");
    }
    input_vector.reserve(255);
    V avg = std::accumulate(t_series_data.begin(), t_series_data.end(), 0,
                            [](auto a, auto b) { return a + b.second; }) /
            (V)t_series_data.size();
    V var = std::accumulate(t_series_data.begin(), t_series_data.end(), 0,
                            [avg](auto a, auto b) {
                                return a + (b.second - avg) * (b.second - avg);
                            }) /
            (V)(t_series_data.size() - 1);
    auto std_dev = std::sqrt(var);
    auto it = t_series_data.begin();
    auto prev = t_series_data.begin();
    input_vector.push_back((prev->second - avg) / std_dev);
    for (it = std::next(it); it != std::end(t_series_data);
         it = std::next(it), prev = std::next(prev)) {
        std::size_t abs_diff = it->first - prev->first;
        if (abs_diff > block_size) {
            for (std::size_t i = 0; i < abs_diff; i += block_size) {
                input_vector.push_back(0);
            }
        }
        input_vector.push_back((prev->second - avg) / std_dev);
    }
    input_vector.shrink_to_fit();

    std::vector<double> in(input_vector);
    // double average = std::accumulate(in.begin(), in.end(), 0) / in.size();
    // for (auto &i : in)
    //   i -= average;

    const std::size_t FFT_ARRAY_SIZE = (in.size() / 2 + 1);

    auto* fft =
        (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * FFT_ARRAY_SIZE);
    fftw_plan forward =
        fftw_plan_dft_r2c_1d(in.size(), in.data(), fft, FFTW_ESTIMATE);
    fftw_execute(forward);
    fftw_destroy_plan(forward);

    double freq_block = 1 / ((double)block_size * in.size());

    std::map<double, std::complex<double>> output_map;
    for (std::size_t i = 0; i < FFT_ARRAY_SIZE; i++) {
        output_map[i * freq_block] = std::complex(fft[i][0], fft[i][1]);
    }
    fftw_free(fft);
    output_map[0] = 0;

    return output_map;
}

}  // namespace fit

template <fit::FitFunction function>
class Chisquared : public ROOT::Minuit2::FCNBase {
private:
    std::vector<double> depths;
    std::vector<double> speed_of_sound;
    std::vector<double> errors;

public:
    Chisquared(std::vector<double> depths, std::vector<double> speed_of_sound,
               std::vector<double> errors)
        : depths{std::move(depths)}
        , speed_of_sound{std::move(speed_of_sound)}
        , errors{std::move(errors)} {};
    ~Chisquared() {};

    double operator()(const std::vector<double>& par) const override;

    [[nodiscard]] std::vector<double>
    fitted_to_minimisation(const ROOT::Minuit2::FunctionMinimum& min) const;
    [[nodiscard]] std::pair<double, double>
    function_minimum(const ROOT::Minuit2::FunctionMinimum& min) const;
    [[nodiscard]] double Up() const override { return 1; };
};

}  // namespace project
