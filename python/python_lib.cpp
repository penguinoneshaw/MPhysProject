#include "speed_of_sound.hpp"
#include <cmath>
#include <vector>

extern "C" {
double unesco_depth(double z, double t, double s, double lat = 0) {
  double pressure = speed_of_sound::pressure_at_depth<double>(z, lat);

  return speed_of_sound::speed_of_sound<double>(pressure, t, s);
}

double unesco_pressure(double pressure, double t, double s) {
  return speed_of_sound::speed_of_sound<double>(pressure, t, s);
}

double leroy_et_al(double z, double t, double s, double lat = 0) {
  return speed_of_sound::leroy_et_al<double>(z, t, s, lat);
}

double ideal_sound_channel(double z, double z_1 = 1160, double B = 1.3,
                           double C_1 = 1.45, double gamma = 1.14e-2) {
  /** https://asa.scitation.org/doi/pdf/10.1121/1.1914492 */
  double eps = 0.5 * B * gamma, eta = (z - z_1) / (1000 * 0.5 * B);
  return C_1 * (1 + eps * (eta + std::exp(-eta) - 1));
}
}