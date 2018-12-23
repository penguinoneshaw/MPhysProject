
#if !defined(H_oQTM)
#define H_oQTM
#include <array>
#include <map>
#include <cmath>
#include <tuple>
#include <memory>

template <typename T, std::size_t LEVEL> class oQTM_Quadrant {
private:
  std::array<std::unique_ptr<oQTM_Quadrant<T, LEVEL + 1> >, 4> subtrees;
  std:map<T, T> data;
public:
  oQTM_Quadrant(T stuff);
  ~oQTM_Quadrant(){};
};

template <typename T, std::size_t N_LEVELS> class oQTM_Mesh {
private:
  std::array<oQTM_Quadrant<T, 1>, 8> quadrants;

public:
  oQTM_Mesh();

  const std::tuple<uint8_t, std::array<uint8_t, N_LEVELS>>
  location(float_t lat, float_t lon) {
    uint8_t octant;
    std::array<uint8_t, N_LEVELS> location;

    if (lat > 0) {
      if (lon < -90) {
        octant = 0;
      } else if (lon < 0) {
        octant = 1;
      } else if (lon < 90) {
        octant = 2;
      } else {
        octant = 3;
      }
    } else {
      if (lon < -90) {
        octant = 4;
      } else if (lon < 0) {
        octant = 5;
      } else if (lon < 90) {
        octant = 6;
      } else {
        octant = 7;
      }
    }

    auto nextlevel = [](const float_t &_x, const float_t &_y) {
      float_t x, y;
      uint8_t quadrant;

      if (_y > 0.5) {
        quadrant = 1;
        x = 2 * _x;
        y = (_y - 0.5) * 2;
      } else if (_y < 0.5 - _x) {
        quadrant = 2;
        x = 2 * _x;
        y = 2 * _y;
      } else if (_x >= 0) {
        y = 2 * _y;
        x = (_x - 0.5) * 2;
        quadrant = 3;
      } else {
        quadrant = 0;
        x = 1 - 2 * _x;
        y = 1 - 2 * _y;
      }

      return std::tie(x, y, quadrant);
    };

    lon = std::fmod(lon + 180.0, 90) / 90.0;
    lat = std::abs(lat) / 90;
    lon *= 1 - lat;

    for (std::size_t i = 0; i < N_LEVELS; i++) {
      auto [nlon, nlat, nlocation] = nextlevel(lon, lat);
      lon = nlon;
      lat = nlat;
      location[i] = nlocation;
    }

    return location;
  }
};

#endif // H_oQTM
