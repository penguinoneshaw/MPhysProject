// #include "./util.hpp"
#include <Minuit2/MnUserParameters.h>

#include <cmath>
#include <vector>

#include "projectfit.hpp"
#include "ranges"

namespace project {

template <>
double
Chisquared<fit::FIT_IDEAL>::operator()(const std::vector<double>& par) const {
    auto fitted_speeds = this->depths | std::views::transform([&par](double z) {
                             return fit::ideal_sound_channel(par, z);
                         });

    double chi2Error{0};
    for (auto [fit, actual, err] :
         std::views::zip(fitted_speeds, this->speed_of_sound, this->errors)) {
        chi2Error += std::pow(fit - actual, 2) / err;
    }
    return chi2Error;
}

template <>
std::pair<double, double> Chisquared<fit::FIT_IDEAL>::function_minimum(
    const ROOT::Minuit2::FunctionMinimum& min) const {
    /**  Calculates the minumum of the function, as well as the error, assuming
     * that a quadratic fitting function was used.
     */
    auto minCoeff = min.UserParameters().Params();

    return {minCoeff[0], min.UserParameters().Error(0)};
}

template <>
std::vector<double> Chisquared<fit::FIT_IDEAL>::fitted_to_minimisation(
    const ROOT::Minuit2::FunctionMinimum& min) const {
    auto par = min.UserState().Params();
    auto output = this->depths | std::views::transform([&par](double d) {
                      return fit::ideal_sound_channel(par, d);
                  });

    return std::ranges::to<std::vector<double>>(output);
}

}  // namespace project

namespace project::fit {
template <>
ROOT::Minuit2::MnUserParameters build_user_parameters<FIT_IDEAL>(
    ROOT::Minuit2::MnUserParameters& user_parameters,
    const double estimated_minimum_depth,
    const double estimated_minimum_speed) {
    user_parameters.Add("z_1", estimated_minimum_depth, 10.0, 0, 100);
    user_parameters.Add("B", 1300, 1.0, 0, 100);
    user_parameters.Add("C_1", estimated_minimum_speed, 1.0, 0, 100);
    user_parameters.Add("gamma", 1.14, 1, 0, 10);

    return user_parameters;
}

double ideal_sound_channel(const std::vector<double>& par, double z) {
    /** https://asa.scitation.org/doi/pdf/10.1121/1.1914492 */
    double z_1 = par[0];
    double B = par[1];
    double C_1 = par[2];
    double gamma = par[3];
    double eps = 0.5 * B * gamma;
    double eta = (z - z_1) / (0.5 * B);
    return C_1 * (1 + eps * (eta + std::expm1(-eta)));
}
}  // namespace project::fit