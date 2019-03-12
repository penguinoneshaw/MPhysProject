#include "speed_of_sound.hpp"

extern "C" {
  double_t unesco_depth(double_t z, double_t t, double_t s, double lat=0){
  double_t pressure = speed_of_sound::pressure_at_depth<double>(z, lat);
    
    return speed_of_sound::speed_of_sound<double> (pressure, t, s);
  };

  double_t unesco_pressure(double_t pressure, double_t t, double_t s){
    return speed_of_sound::speed_of_sound<double> (pressure, t, s);

  }

  double_t leroy_et_al(double_t z, double_t t, double_t s, double lat=0){
    return speed_of_sound::leroy_et_al<double>(z,t,s,lat);
  }

}