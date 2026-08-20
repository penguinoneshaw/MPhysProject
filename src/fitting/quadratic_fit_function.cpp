// #include "./util.hpp"
#include <Minuit2/MnUserParameters.h>

#include <cmath>
#include <vector>

#include "polynomial.hpp"
#include "projectfit.hpp"
#include "ranges"

namespace project {

double polynomial_fit(const std::vector<double>& par, double d) {
    return poly::horners_method(par, d);
}

template <>
double Chisquared<fit::FIT_QUADRATIC>::operator()(
    const std::vector<double>& par) const {
    auto fitted_speeds =
        std::views::zip(this->depths | std::views::transform([&par](double d) {
                            return polynomial_fit(par, d);
                        }),
                        this->speed_of_sound, this->errors);

    double result = 0;

    for (auto [fit, actual, err] : fitted_speeds) {
        result += std::pow(fit - actual, 2) / err;
    }
    return result;
}

template <>
std::pair<double, double> Chisquared<fit::FIT_QUADRATIC>::function_minimum(
    const ROOT::Minuit2::FunctionMinimum& min) const {
    /**  Calculates the minumum of the function, as well as the error, assuming
     * that a quadratic fitting function was used.
     */
    auto minCoeff = min.UserParameters().Params();
    return {
        -minCoeff[1] / (2 * minCoeff[2]),
        std::abs(-minCoeff[1] / (2 * minCoeff[2])) *
            std::sqrt(std::pow(min.UserParameters().Error(1) / minCoeff[1], 2) +
                      std::pow(min.UserParameters().Error(2) / minCoeff[2], 2))};
}

template <>
std::vector<double> Chisquared<fit::FIT_QUADRATIC>::fitted_to_minimisation(
    const ROOT::Minuit2::FunctionMinimum& min) const {
    auto par = min.UserState().Params();
    auto output = this->depths | std::views::transform([&par](double d) {
                      return polynomial_fit(par, d);
                  });

    return std::ranges::to<std::vector<double>>(output);
}

}  // namespace project

namespace project::fit {
template <>
ROOT::Minuit2::MnUserParameters build_user_parameters<FIT_QUADRATIC>(
    ROOT::Minuit2::MnUserParameters& user_parameters,
    const double estimated_minimum_depth,
    const double estimated_minimum_speed) {
    user_parameters.Add("c_0", estimated_minimum_depth, 1);
    user_parameters.Add("c_1", 0, 10);
    user_parameters.Add("c_2", 0, 10);

    return user_parameters;
}
}  // namespace project::fit