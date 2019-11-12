// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2019 Scipp contributors (https://github.com/scipp)
/// @file
/// @author Simon Heybrock
#ifndef SCIPP_CORE_TRANSFORM_COMMON_H
#define SCIPP_CORE_TRANSFORM_COMMON_H

#include <tuple>

namespace scipp::core {

template <class... Ts> struct pair_self {
  using type = std::tuple<std::pair<Ts, Ts>...>;
};
template <class... Ts> struct pair_custom { using type = std::tuple<Ts...>; };
template <class... Ts> struct pair_ {
  template <class RHS> using type = std::tuple<std::pair<Ts, RHS>...>;
};

template <class... Ts> using pair_self_t = typename pair_self<Ts...>::type;
template <class... Ts> using pair_custom_t = typename pair_custom<Ts...>::type;
template <class RHS>
using pair_numerical_with_t =
    typename pair_<double, float, int64_t, int32_t>::type<RHS>;

template <class... Ts> struct pair_product {
  template <class T> using type = std::tuple<std::pair<T, Ts>...>;
};

template <class... Ts>
using pair_product_t = decltype(
    std::tuple_cat(typename pair_product<Ts...>::template type<Ts>{}...));

using arithmetic_type_pairs = pair_product_t<float, double, int32_t, int64_t>;

using arithmetic_type_pairs_with_bool =
    decltype(std::tuple_cat(std::declval<arithmetic_type_pairs>(),
                            std::declval<pair_numerical_with_t<bool>>()));

using arithmetic_and_matrix_type_pairs = decltype(std::tuple_cat(
    std::declval<arithmetic_type_pairs>(),
    std::tuple<std::pair<Eigen::Vector3d, Eigen::Vector3d>,
               std::pair<int64_t, int32_t>, std::pair<int32_t, int64_t>,
               std::pair<double, float>, std::pair<float, double>>()));

static constexpr auto dimensionless_unit_check =
    [](units::Unit &varUnit, const units::Unit &otherUnit) {
      expect::equals(varUnit, units::dimensionless);
      expect::equals(otherUnit, units::dimensionless);
    };

static constexpr auto dimensionless_unit_check_return =
    [](const units::Unit &aUnit, const units::Unit &bUnit) {
      expect::equals(aUnit, units::dimensionless);
      expect::equals(bUnit, units::dimensionless);
      return aUnit;
    };

} // namespace scipp::core

#endif // SCIPP_CORE_TRANSFORM_COMMON_H