// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2020 Scipp contributors (https://github.com/scipp)
/// @file
/// @author Owen Arnold
#ifndef SCIPP_CORE_ELEMENT_GEOMETRIC_OPERATIONS_H
#define SCIPP_CORE_ELEMENT_GEOMETRIC_OPERATIONS_H

#include "arg_list.h"
#include "scipp/common/overloaded.h"
#include "scipp/core/transform.h"

namespace scipp::core {

/// Operators to be used with transform and transform_in_place to implement
/// geometric operations for Variable.

namespace element {

namespace geometry {

constexpr auto position = overloaded{
    arg_list<double>,
    transform_flags::expect_no_variance_arg<0>,
    transform_flags::expect_no_variance_arg<1>,
    transform_flags::expect_no_variance_arg<2>,
    [](const auto &x, const auto &y, const auto &z) {
      using T = double; // currently only double precision support
      return Eigen::Matrix<T, 3, 1>(x, y, z);
    },
    [](const units::Unit &x, const units::Unit &y, const units::Unit &z) {
      expect::equals(x, y);
      expect::equals(x, z);
      expect::equals(x, units::m);
      return x;
    }};

namespace detail {
template <int N> struct component {
  static constexpr auto overloads = overloaded{
      arg_list<Eigen::Vector3d>, [](const auto &pos) { return pos[N]; },
      [](const units::Unit &u) {
        expect::equals(u, units::m);
        return u;
      }};
  enum { value = N };
};
} // namespace detail
constexpr auto x = detail::component<0>::overloads;
constexpr auto y = detail::component<1>::overloads;
constexpr auto z = detail::component<2>::overloads;

constexpr auto rotate =
    overloaded{
      arg_list<std::tuple<Eigen::Vector3d, Eigen::Quaterniond>>,
      transform_flags::expect_no_variance_arg<0>,
      transform_flags::expect_no_variance_arg<1>,
      [](const auto &pos, const auto &rot) {
                 return rot._transformVector(pos);
               },
      [](const units::Unit &u_pos, const units::Unit &u_rot) {
        expect::equals(u_pos, units::m);
        expect::equals(u_rot, units::dimensionless);
        return u_pos;
      }};
} // namespace geometry

} // namespace element

} // namespace scipp::core

#endif // SCIPP_CORE_ELEMENT_GEOMETRIC_OPERATIONS_H
