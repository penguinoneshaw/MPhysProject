#include "projectfit.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <memory>
#include <numeric>
#include <ranges>
#include <vector>

#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnUserParameters.h"
#include "fftw3.h"
#include "util.hpp"

namespace project {

double bilinear_fit(const std::vector<double>& par, double d) {
    // This is the fitting function, which I can swap out as necessary
    // Uses a bilinear fit, where the curve is theoretically fitted to two
    // successive lines, which can be generalised. Expects par to be in the format
    // (x, y, m_1, m_2) where [x, y] is the 'position' of the turning point in the
    // depth-speed plane

    double r = 0;
    if (d < par[0]) {
        r = par[2] * (d - par[0]) + par[1];
    } else {
        r = par[3] * (d - par[0]) + par[1];
    }
    return r;
}

namespace fit {
std::vector<float> low_pass_filter(const std::vector<float_t>& vector,
                                   const std::size_t cutoff) {
    float average = project::utils::average(vector);

    auto in = std::ranges::to<std::vector<float_t>>(
        vector |
        std::views::transform([average](float_t v) { return v - average; }));

    const std::size_t FFT_ARRAY_SIZE = (in.size() / 2 + 1);
    std::vector<float> out(in.size());
    fftwf_complex* fft =
        (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * FFT_ARRAY_SIZE);
    fftwf_plan forward =
        fftwf_plan_dft_r2c_1d(in.size(), in.data(), fft, FFTW_ESTIMATE);
    fftwf_execute(forward);
    fftwf_destroy_plan(forward);

    for (std::size_t i = 0; i < FFT_ARRAY_SIZE; i++) {
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

    for (std::size_t i = 0; i < out.size(); i++) {
        out[i] /= out.size();
        out[i] += average;
    }

    return out;
}

std::vector<double> low_pass_filter(const std::vector<double>& vector,
                                    const std::size_t cutoff) {
    float average = project::utils::average(vector);

    auto in = std::ranges::to<std::vector<double>>(
        vector |
        std::views::transform([average](float v) { return v - average; }));

    std::vector<double> out(in.size());
    fftw_complex* fft =
        (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * in.size() / 2 + 1);
    fftw_plan forward =
        fftw_plan_dft_r2c_1d(in.size(), in.data(), fft, FFTW_ESTIMATE);
    fftw_execute(forward);
    fftw_destroy_plan(forward);

    for (std::size_t i = 0; i < in.size() / 2 + 1; i++) {
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

    for (std::size_t i = 0; i < out.size(); i++) {
        out[i] /= out.size();
        out[i] += average;
    }

    return out;
}

template <std::floating_point T>
std::vector<T> low_pass_filter(const std::vector<T>& vector,
                               const std::size_t cutoff) {
    if constexpr (std::is_same_v<T, float_t>) {
        return low_pass_filter(std::ranges::to<std::vector<float_t>>(vector),
                               cutoff);
    } else if constexpr (std::is_same_v<T, double>) {
        return low_pass_filter(std::ranges::to<std::vector<double>>(vector),
                               cutoff);
    } else {
        throw std::runtime_error("Type not supported for low pass filter");
    }
}

template <std::floating_point T>
std::vector<T> differentiate(const std::vector<T>& xs,
                             const std::vector<T>& ys) {
    std::vector<T> result(ys.size(), 0), dxs(xs.size(), 0);
    std::adjacent_difference(ys.begin(), ys.end(), result.begin());
    std::adjacent_difference(xs.begin(), xs.end(), dxs.begin());
    for (size_t i = 0; i < result.size() && i < dxs.size(); i++) {
        result[i] = result[i] / dxs[i];
    }

    return result;
}

template <std::ranges::input_range T>
    requires std::floating_point<std::ranges::range_value_t<T>>
std::ranges::range_value_t<T> global_minimum(const T& range) {
    return *std::ranges::min_element(range);
}

template <std::floating_point T, FitFunction f = FitFunction::FIT_QUADRATIC>
std::pair<T, T> find_sofar_channel(const std::vector<T>& speed_of_sound,
                                   const std::vector<T>& depths,
                                   std::size_t averaging_granularity = 10) {
    auto [avg_sos, sos_errors] =
        project::utils::moving_average(speed_of_sound, averaging_granularity);

    auto [avg_depths, depths_errors] =
        project::utils::moving_average(depths, averaging_granularity);

    auto diff_avg_sos = differentiate(avg_depths, avg_sos);

    std::vector<size_t> maxima{0};
    if (diff_avg_sos.size() < 5) {
        throw std::runtime_error("NOT ENOUGH DATA");
    }

    for (auto it = diff_avg_sos.begin() + 1; it != diff_avg_sos.end(); ++it) {
        if (*(it - 1) > 0 && *it <= 0) {
            maxima.push_back(std::distance(diff_avg_sos.begin(), it));
        }
    }

    std::size_t threshold = avg_sos.size() > 10 ? avg_sos.size() - 10 : 0;
    std::size_t minimum_index = maxima.back() >= threshold ? 0 : maxima.back();

    std::vector<double> chisquared_depths(std::begin(avg_depths) + minimum_index,
                                          avg_depths.end());
    auto chisquared_sos = avg_sos | std::views::drop(minimum_index);

    auto chisquared_errors = sos_errors | std::views::drop(minimum_index);

    Chisquared<f> fcn(chisquared_depths,
                      std::ranges::to<std::vector<double>>(chisquared_sos),
                      std::ranges::to<std::vector<double>>(chisquared_errors));

    ROOT::Minuit2::MnUserParameters user_parameters{};
    build_user_parameters<f>(user_parameters, depths[minimum_index],
                             speed_of_sound[minimum_index]);

    ROOT::Minuit2::MnMigrad migrad(fcn, user_parameters);
    ROOT::Minuit2::FunctionMinimum min = migrad();
    auto [xmin, err] = fcn.function_minimum(min);

    if (!min.IsValid() || std::isnan(err) || std::isnan(xmin) ||
        xmin > avg_depths.back() || xmin < avg_depths.front()
        /*|| err > 20*/) {
        throw std::runtime_error("Fitted minimum invalid, probably out of region");
    }

    return std::make_pair(static_cast<T>(depths[minimum_index]),
                          static_cast<T>(depths_errors[minimum_index]));
}

}  // namespace fit

}  // end namespace project
template std::pair<double, double>
project::fit::find_sofar_channel<double,
                                 project::fit::FitFunction::FIT_QUADRATIC>(
    const std::vector<double>&, const std::vector<double>&, std::size_t);
template std::pair<double, double>
project::fit::find_sofar_channel<double, project::fit::FitFunction::FIT_IDEAL>(
    const std::vector<double>&, const std::vector<double>&, std::size_t);