
#if !defined(H_oQTM)
#define H_oQTM

#include <array>
#include <map>
#include <cmath>
#include <tuple>
#include <memory>

template <typename T, std::size_t N_LEVELS>
class oQTM_Quadrant
{
private:
  typedef oQTM_Quadrant<T, N_LEVELS> subtree;
  typedef std::unique_ptr<subtree> subtree_ptr;

  class BeyondTreeDepthError : public std::exception
  {
  public:
    virtual char const *what() const noexcept { return "Attempting to access level of tree beyond that in template type."; }
  } beyondTreeError;

  std::array<subtree_ptr, 4> subtrees;
  std::map<T, T> data;
  std::size_t depth;

public:
  oQTM_Quadrant(std::size_t _depth = 1): depth{_depth} {

  };
  ~oQTM_Quadrant(){};
  subtree &operator[](std::size_t i)
  {
    if (depth == N_LEVELS) {
      throw beyondTreeError;
    }
    if (!subtrees[i])
      subtrees[i] = std::make_unique<subtree>(depth + 1);
    return *subtrees[i];
  }
};

template <typename T, std::size_t N_LEVELS>
class oQTM_Mesh
{
private:
  typedef oQTM_Quadrant<T, N_LEVELS> octant;
  typedef std::unique_ptr<octant> quad_ptr;
  std::array<quad_ptr, 8> quadrants;

public:
  oQTM_Mesh(){};
  ~oQTM_Mesh(){};

  octant &operator[](std::size_t i)
  {
    if (!quadrants[i])
      quadrants[i] = std::make_unique<octant>();
    return *quadrants[i];
  }

  static std::tuple<uint8_t, uint8_t, T> nextlevel(const T &_x, const T &_y)
  {
    T x, y;
    uint8_t quadrant;

    if (_y > 0.5)
    {
      quadrant = 1;
      x = 2 * _x;
      y = (_y - 0.5) * 2;
    }
    else if (_y < 0.5 - _x)
    {
      quadrant = 2;
      x = 2 * _x;
      y = 2 * _y;
    }
    else if (_x >= 0)
    {
      y = 2 * _y;
      x = (_x - 0.5) * 2;
      quadrant = 3;
    }
    else
    {
      quadrant = 0;
      x = 1 - 2 * _x;
      y = 1 - 2 * _y;
    }

    return std::tie(x, y, quadrant);
  };

  const std::tuple<uint8_t, std::array<uint8_t, N_LEVELS>>
  location(float_t lat, float_t lon) const
  {
    uint8_t octant;
    std::array<uint8_t, N_LEVELS> location;

    if (lat > 0)
    {
      if (lon < -90)
      {
        octant = 0;
      }
      else if (lon < 0)
      {
        octant = 1;
      }
      else if (lon < 90)
      {
        octant = 2;
      }
      else
      {
        octant = 3;
      }
    }
    else
    {
      if (lon < -90)
      {
        octant = 4;
      }
      else if (lon < 0)
      {
        octant = 5;
      }
      else if (lon < 90)
      {
        octant = 6;
      }
      else
      {
        octant = 7;
      }
    }

    lon = std::fmod(lon + 180.0, 90) / 90.0;
    lat = std::abs(lat) / 90;
    lon *= 1 - lat;

    for (std::size_t i = 0; i < N_LEVELS; i++)
    {
      auto [nlon, nlat, nlocation] = nextlevel(lon, lat);
      lon = nlon;
      lat = nlat;
      location[i] = nlocation;
    }

    return location;
  }
};

#endif // H_oQTM
