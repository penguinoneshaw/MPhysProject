
#if !defined(H_oQTM)
#define H_oQTM

#include <array>
#include <map>
#include <cmath>
#include <tuple>
#include <memory>
#include <stdexcept>

class BeyondTreeDepthError : public std::logic_error
{
public:
  BeyondTreeDepthError() : std::logic_error("Attempting to access level of tree beyond that in template type.")
  {
  }
} beyondTreeError;

template <typename K, typename V, std::size_t N_LEVELS>
class oQTM_Quadrant
{
private:
  typedef oQTM_Quadrant<K, V, N_LEVELS> subtree;
  typedef std::shared_ptr<subtree> subtree_ptr;

  std::array<subtree_ptr, 4> subtrees;
  std::multimap<K, V> data;
  std::size_t depth;

public:
  oQTM_Quadrant(std::size_t _depth = 1) : depth{_depth} {};
  ~oQTM_Quadrant(){};

  void add_to_data(K key, V value) {
    data.insert(std::pair(key, value));
  }

  const std::multimap<K, V> get_data()
  {
    if (depth == N_LEVELS)
    {
      return data;
    }
    else
    {
      std::multimap<K, V> combined_map;
      for (auto &subtree : subtrees)
      {
        if (subtree)
        {
          auto st = subtree->get_data();
          combined_map.insert(st.begin(), st.end());
        }
      }
      return combined_map;
    }
  }

  std::shared_ptr<subtree> operator[](std::size_t i)
  {
    if (depth >= N_LEVELS)
    {
      throw beyondTreeError;
    }

    if (!subtrees[i])
      subtrees[i] = std::make_shared<subtree>(depth + 1);
    return subtrees[i];
  }

  std::multimap<K, V> get_points(const std::vector<size_t>::iterator begin, const std::vector<size_t>::const_iterator &end)
  {
    if (begin == end)
    {
      return this->get_data();
    }
    else
    {
      return this->operator[](*begin)->get_points(begin + 1, end);
    }
  }
};

template <typename T, typename K, typename V, std::size_t N_LEVELS>
class oQTM_Mesh
{
  private:
  typedef oQTM_Quadrant<K, V, N_LEVELS> octant;
  typedef std::shared_ptr<octant> quad_ptr;
  std::array<quad_ptr, 8> quadrants;

  public:
  oQTM_Mesh(){};
  ~oQTM_Mesh(){};
  typedef std::array<uint8_t, N_LEVELS> location_t;

  std::shared_ptr<octant> operator[](std::size_t i)
  {
    if (!quadrants[i])
      quadrants[i] = std::make_shared<octant>();
    return quadrants[i];
  }




  std::multimap<K, V> get_points(std::vector<size_t> location)
  {
    if (location.size() > N_LEVELS)
      throw beyondTreeError;
    auto begin = location.begin();
    return operator[](*begin)->get_points(begin + 1, location.end());
  }

  static inline std::tuple<T, T, uint8_t> nextlevel(const T &_x, const T &_y)
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

  const location_t
  location(T lon, T lat) const
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

    for (std::size_t i = 1; i < N_LEVELS; i++)
    {
      auto [nlon, nlat, nlocation] = nextlevel(lon, lat);
      lon = nlon;
      lat = nlat;
      location[i] = nlocation;
    }

    location[0] = octant;

    return location;
  }

  const location_t insert(T lon, T lat, K key, V value){
    auto loc = location(lon, lat);
    this->insert(loc, key, value);
    return loc;
  }

  void insert(location_t location, K key, V value){
    std::shared_ptr<octant> quad = operator[](location[0]);
    for (auto it = location.begin() + 1; it != location.end(); it++){
      quad = (*quad)[*it];
    }
    quad->add_to_data(key, value);
  }
};

#endif // H_oQTM
