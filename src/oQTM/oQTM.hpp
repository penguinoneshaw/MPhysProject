#pragma once

#include <array>
#include <cassert>
#include <cmath>
#include <concepts>
#include <map>
#include <memory>
#include <mutex>
#include <ranges>
#include <stdexcept>
#include <tuple>
#include <vector>
class BeyondTreeDepthError : public std::logic_error {
public:
    BeyondTreeDepthError()
        : std::logic_error(
              "Attempting to access level of tree beyond that in "
              "template type.") {}
};
const BeyondTreeDepthError beyondTreeError{};

template <typename K, typename V, std::size_t N_LEVELS>
class oQTM_Quadrant {
private:
    using subtree = oQTM_Quadrant<K, V, N_LEVELS>;
    using subtree_ptr = std::shared_ptr<subtree>;
    using location_t = std::array<uint8_t, N_LEVELS>;
    std::array<subtree_ptr, 4> subtrees;
    std::multimap<K, V> data;
    std::size_t depth{};

public:
    explicit oQTM_Quadrant(std::size_t _depth = 1)
        : depth{_depth} {};

    oQTM_Quadrant(const subtree& quadrant) {
        if (quadrant.is_datalayer()) {
            this->data = quadrant.get_data();
        } else {
            auto quad_subtrees = quadrant.get_subtrees();
            for (std::size_t i{0}; i < 4; i++) {
                if (quad_subtrees[i]) {
                    this->subtrees[i] = std::make_shared<oQTM_Quadrant>(*quad_subtrees[i]);
                }
            }
        }
    };
    ~oQTM_Quadrant() = default;

    std::array<subtree_ptr, 4> get_subtrees() const { return subtrees; }

    [[nodiscard]] bool is_datalayer() const { return depth == N_LEVELS; }

    void add_to_data(K key, V value) { data.insert(std::pair(key, value)); }

    std::multimap<K, V> get_data() const {
        if (is_datalayer()) {
            return data;
        } else {
            std::multimap<K, V> combined_map;
            for (auto& subtree : subtrees) {
                if (subtree) {
                    auto st = subtree->get_data();
                    combined_map.insert(st.begin(), st.end());
                }
            }
            return combined_map;
        }
    }

    std::shared_ptr<subtree> operator[](std::size_t i) {
        if (depth >= N_LEVELS) {
            throw BeyondTreeDepthError();
        }

        if (!subtrees[i])
            subtrees[i] = std::make_shared<subtree>(depth + 1);
        return subtrees[i];
    }

    std::multimap<K, V>
    get_points(const typename location_t::iterator begin,
               const typename location_t::const_iterator& end) {
        if (begin == end) {
            return this->get_data();
        } else {
            return this->operator[](*begin)->get_points(begin + 1, end);
        }
    }
};

template <std::floating_point TCoordinate, typename K, typename V,
          std::size_t N_LEVELS, bool ERROR = false>
class oQTM_Mesh {
private:
    using octant = oQTM_Quadrant<K, V, N_LEVELS>;
    using quad_ptr = std::shared_ptr<octant>;
    std::array<quad_ptr, 8> octants;
    std::mutex write_mutex;

public:
    oQTM_Mesh() { octants.fill(nullptr); };

    explicit oQTM_Mesh(const oQTM_Mesh<TCoordinate, K, V, N_LEVELS>& mesh) {
        for (std::size_t i = 0; i < 8; i++) {
            if (mesh.get_subtree_without_creating(i)) {
                this->octants[i] = std::make_shared<octant>(*mesh.get_subtree_without_creating(i));
            }
        }
    }

    ~oQTM_Mesh() = default;

    using location_t = std::array<uint8_t, N_LEVELS>;

    std::shared_ptr<octant> operator[](std::size_t i) {
        if (!octants[i])
            octants[i] = std::make_shared<octant>();
        return octants[i];
    }

    inline const quad_ptr& get_subtree_without_creating(std::size_t i) const {
        return octants[i];
    }

    std::multimap<K, V> get_points(location_t location,
                                   std::size_t granularity = N_LEVELS) {
        if (location.size() > N_LEVELS) {
            throw BeyondTreeDepthError();
        }
        auto begin = location.begin();
        return operator[](*begin)->get_points(
            begin + 1, granularity < N_LEVELS ? location.begin() + granularity
                                              : location.end());
    }

    std::map<K, V> get_averaged_points(location_t location,
                                       std::size_t spatial_granularity = N_LEVELS,
                                       std::size_t temporal_granularity = 1) {
        assert(temporal_granularity != 0);
        if (location.size() > N_LEVELS)
            throw BeyondTreeDepthError();

        std::map<K, V> result;
        std::map<K, uint64_t> count;
        for (auto& element : get_points(location, spatial_granularity)) {
            K index = std::floor(std::round(element.first / temporal_granularity) *
                                 temporal_granularity);
            ++count[index];
            if (ERROR) {
                result[index] += std::pow(element.second, 2);
            } else {
                result[index] += element.second;
            }
        }

        for (auto& [element, _] : result) {
            result[element] /= ((V)count[element]);
            if (ERROR) {
                result[element] = std::sqrt(result[element]);
            }
        }

        return result;
    }

    static inline std::tuple<TCoordinate, TCoordinate, uint8_t>
    nextlevel(const TCoordinate& _x, const TCoordinate& _y) {
        TCoordinate x;
        TCoordinate y;
        uint8_t quadrant;

        if (_y > 0.5) {
            quadrant = 1;
            x = 2.0 * _x;
            y = (_y - 0.5) * 2;
        } else if (_y < 0.5 - _x) {
            quadrant = 2;
            x = 2.0 * _x;
            y = 2.0 * _y;
        } else if (_x >= 0) {
            y = 2.0 * _y;
            x = (_x - 0.5) * 2.0;
            quadrant = 3;
        } else {
            quadrant = 0;
            x = 1.0 - 2.0 * _x;
            y = 1.0 - 2.0 * _y;
        }

        return std::make_tuple(x, y, quadrant);
    };

    constexpr location_t location(TCoordinate lon, TCoordinate lat) const {
        uint8_t octant;
        std::array<uint8_t, N_LEVELS> location;

        if (lat > 0.0) {
            if (lon < -90.0) {
                octant = 0;
            } else if (lon < 0.0) {
                octant = 1;
            } else if (lon < 90.0) {
                octant = 2;
            } else {
                octant = 3;
            }
        } else {
            if (lon < -90.0) {
                octant = 4;
            } else if (lon < 0.0) {
                octant = 5;
            } else if (lon < 90.0) {
                octant = 6;
            } else {
                octant = 7;
            }
        }

        lon = std::fmod(lon + 180.0, 90.0) / 90.0;
        lat = std::abs(lat) / 90.0;
        lon *= 1.0 - lat;

        for (std::size_t i = 1; i < N_LEVELS; i++) {
            auto [nlon, nlat, nlocation] = nextlevel(lon, lat);
            lon = nlon;
            lat = nlat;
            location[i] = nlocation;
        }

        location[0] = octant;

        return location;
    }

    location_t insert(TCoordinate lon, TCoordinate lat, K key, V value) {
        auto loc = location(lon, lat);
        this->insert(loc, key, value);
        return loc;
    }

    void insert(location_t location, K key, V value) {
        std::shared_ptr<octant> quad = operator[](location[0]);
        for (auto locationIndex : std::views::drop(location, 1)) {
            quad = quad->operator[](locationIndex);
        }
        quad->add_to_data(key, value);
    }

    void insert(
        const std::vector<std::tuple<TCoordinate, TCoordinate, K, V>> locations) {
        std::lock_guard lock(write_mutex);
        for (auto point : locations) {
            auto [lon, lat, key, value] = point;
            this->insert(lon, lat, key, value);
        }
    }
};
