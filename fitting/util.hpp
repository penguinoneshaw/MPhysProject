#pragma once
#include "iterator"
#include <numeric>
#include <ranges>
#include <tuple>
#include <utility>
#include <vector>

template <std::ranges::input_range... Ranges>
  requires(sizeof...(Ranges) > 0)
class zip_view : public std::ranges::view_interface<zip_view<Ranges...>> {
private:
  std::tuple<Ranges...> bases;

  template <typename... Its> struct zip_iterator {
    std::tuple<Its...> its;

    using value_type = std::tuple<std::ranges::range_value_t<Ranges>...>;
    using reference = std::tuple<std::ranges::range_reference_t<Ranges>...>;
    using iterator_category = std::input_iterator_tag;

    reference operator*() const {
      return std::apply(
          [](auto &...it) { return std::tuple<decltype(*it)...>(*it...); },
          its);
    }

    zip_iterator &operator++() {
      std::apply([](auto &...it) { (++it, ...); }, its);
      return *this;
    }

    bool operator==(const zip_iterator &other) const {
      return std::apply(
          [&](auto &...it) {
            return ((it == std::get<&it - &std::get<0>(its)>(other.its)) ||
                    ...);
          },
          its);
    }
  };

public:
  zip_view(Ranges... rs) : bases(std::move(rs)...) {}

  auto begin() {
    return std::apply(
        [](auto &...r) {
          return zip_iterator<decltype(std::ranges::begin(r))...>{
              std::make_tuple(std::ranges::begin(r)...)};
        },
        bases);
  }

  auto end() {
    return std::apply(
        [](auto &...r) {
          return zip_iterator<decltype(std::ranges::end(r))...>{
              std::make_tuple(std::ranges::end(r)...)};
        },
        bases);
  }
};

namespace project::utils {

template <std::ranges::input_range Range>
  requires std::floating_point<std::ranges::range_value_t<Range>>
auto average(const Range &range) {
  using T = std::ranges::range_value_t<Range>;
  auto first = std::ranges::begin(range);
  auto last = std::ranges::end(range);
  auto count = std::ranges::distance(first, last);
  return std::accumulate(first, last, T{}) / static_cast<T>(count);
}

template <std::ranges::input_range Range>
  requires std::floating_point<std::ranges::range_value_t<Range>>
auto variance(const Range &range) {
  using T = std::ranges::range_value_t<Range>;
  auto first = std::ranges::begin(range);
  auto last = std::ranges::end(range);
  auto count = std::ranges::distance(first, last);

  auto mean = std::accumulate(first, last, T{}) / static_cast<T>(count);
  T sq_sum = std::accumulate(first, last, T{}, [mean](T acc, T val) {
    return acc + (val - mean) * (val - mean);
  });
  return sq_sum / static_cast<T>(count);
}

template <std::ranges::input_range T>
  requires std::floating_point<std::ranges::range_value_t<T>>
std::pair<std::vector<std::ranges::range_value_t<T>>,
          std::vector<std::ranges::range_value_t<T>>>
moving_average(const T &vector, const std::size_t period) {
  auto means = vector | std::views::slide(period) |
               std::views::transform([](const auto &r) { return average(r); });

  auto errors =
      vector | std::views::slide(period) |
      std::views::transform([](const auto &r) { return variance(r); });

  return {std::ranges::to<std::vector<std::ranges::range_value_t<T>>>(means),
          std::ranges::to<std::vector<std::ranges::range_value_t<T>>>(errors)};
}

} // namespace project::utils