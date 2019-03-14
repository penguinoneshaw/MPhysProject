#include "speed_of_sound.hpp"
#include <vector>

extern "C" {
  double_t unesco_depth(double_t z, double_t t, double_t s, double lat=0){
  double_t pressure = speed_of_sound::pressure_at_depth<double>(z, lat);
    
    return speed_of_sound::speed_of_sound<double> (pressure, t, s);
  }

  double_t unesco_pressure(double_t pressure, double_t t, double_t s){
    return speed_of_sound::speed_of_sound<double> (pressure, t, s);

  }

  double_t leroy_et_al(double_t z, double_t t, double_t s, double lat=0){
    return speed_of_sound::leroy_et_al<double>(z,t,s,lat);
  }

  double_t ideal_sound_channel( double_t z, double_t z_1 = 1160, double_t B=1.3, double_t C_1=1.45, double_t gamma=1.14e-2){
  /** https://asa.scitation.org/doi/pdf/10.1121/1.1914492 */
    double eps = 0.5*B*gamma, eta = (z - z_1)/(1000*0.5*B);
    return C_1 * (1 + eps * (eta + std::exp(-eta) - 1));
  }

}