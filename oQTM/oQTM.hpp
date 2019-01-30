
#if !defined(H_oQTM)
#define H_oQTM

#include <array>
#include <map>
#include <cmath>
#include <tuple>
#include <memory>
#include <stdexcept>
#include <mutex>

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
  typedef std::array<uint8_t, N_LEVELS> location_t;
  std::array<subtree_ptr, 4> subtrees;
  std::multimap<K, V> data;
  std::size_t depth;

public:
  oQTM_Quadrant(std::size_t _depth = 1) : depth{_depth} {};

  oQTM_Quadrant(const subtree &quadrant)
  {
    if (quadrant.is_datalayer())
    {
      this->data = quadrant.get_data();
    }
    else
    {
      auto quad_subtrees = quadrant.get_subtrees();
      for (std::size_t i; i < 4; i++)
      {
        if (quad_subtrees[i])
        {
          this->subtrees[i] = oQTM_Quadrant(*quad_subtrees[i]);
        }
      }
    }
  };
  ~oQTM_Quadrant(){};

  const std::array<subtree_ptr, 4> get_subtrees()
  {
    return subtrees;
  }

  bool is_datalayer() const
  {
    return depth == N_LEVELS;
  }
  void add_to_data(K key, V value)
  {
    data.insert(std::pair(key, value));
  }

  const std::multimap<K, V> get_data()
  {
    if (is_datalayer())
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

  std::multimap<K, V> get_points(const typename location_t::iterator begin, const typename location_t::const_iterator &end)
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
  std::array<quad_ptr, 8> octants;
  std::mutex write_mutex;

public:
  oQTM_Mesh(){
    octants.fill(nullptr);
  };

  oQTM_Mesh(const oQTM_Mesh<T, K, V, N_LEVELS> &mesh)
  {
    for (std::size_t i = 0; i < 8; i++)
    {
      if (mesh.get_subtree_without_creating(i))
      {
        this->octants[i] = octant(*mesh.get_subtree_without_creating(i));
      }
    }
  }

  ~oQTM_Mesh(){};
  typedef std::array<uint8_t, N_LEVELS> location_t;

  std::shared_ptr<octant> operator[](std::size_t i)
  {
    if (!octants[i])
      octants[i] = std::make_shared<octant>();
    return octants[i];
  }

  inline const quad_ptr &get_subtree_without_creating(std::size_t i)
  {
    return octants[i];
  }

  std::multimap<K, V> get_points(location_t location, std::size_t granularity = N_LEVELS)
  {
    if (location.size() > N_LEVELS)
      throw beyondTreeError;
    auto begin = location.begin();
    return operator[](*begin)->get_points(begin + 1, granularity < N_LEVELS ? location.begin() + granularity : location.end());
  }

  std::map<K, V> get_averaged_points(location_t location, std::size_t granularity = N_LEVELS)
  {
    if (location.size() > N_LEVELS)
      throw beyondTreeError;

    std::map<K, V> result;
    std::map<K, uint64_t> count;
    for (auto &element : get_points(location, granularity))
    {
      try
      {
        auto curr = result.at(element.first);
        auto curr_count = ++count[element.first];
        result[element.first] = (curr * ((V) (curr_count - 1)) + element.second) / ((V)curr_count);
      }
      catch (std::out_of_range err)
      {
        result[element.first] = element.second;
      }
    }

    return result;
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

    if (lat > 0.0)
    {
      if (lon < -90.0)
      {
        octant = 0;
      }
      else if (lon < 0.0)
      {
        octant = 1;
      }
      else if (lon < 90.0)
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
      if (lon < -90.0)
      {
        octant = 4;
      }
      else if (lon < 0.0)
      {
        octant = 5;
      }
      else if (lon < 90.0)
      {
        octant = 6;
      }
      else
      {
        octant = 7;
      }
    }

    lon = std::fmod(lon + 180.0, 90.0) / 90.0;
    lat = std::abs(lat) / 90.0;
    lon *= 1.0 - lat;

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

  const location_t insert(T lon, T lat, K key, V value)
  {
    auto loc = location(lon, lat);
    this->insert(loc, key, value);
    return loc;
  }

  void insert(location_t location, K key, V value)
  {
    std::shared_ptr<octant> quad = operator[](location[0]);
    for (auto it = location.begin() + 1; it != location.end(); it++)
    {
      quad = (*quad)[*it];
    }
    quad->add_to_data(key, value);
  }

  void insert(const std::vector<std::tuple<T, T, K, V>> locations)
  {
    std::lock_guard lock(write_mutex);
    std::for_each(locations.begin(), locations.end(), [this](const std::tuple<T, T, K, V> &point) -> void {
      auto [lon, lat, key, value] = point;
      this->insert(lon, lat, key, value);
    });
  }
};

#endif // H_oQTM
