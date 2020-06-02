// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2020 Scipp contributors (https://github.com/scipp)
/// @file
/// @author Simon Heybrock
#pragma once

#include "scipp/common/overloaded.h"
#include "scipp/common/span.h"
#include "scipp/core/element/arg_list.h"
#include "scipp/core/value_and_variance.h"
#include "scipp/units/except.h"
#include "scipp/units/unit.h"

#include <stddef.h>
#include <limits>

namespace scipp::core::element {

/// Sets any masked elements to 0 to handle special FP vals
constexpr auto convertMaskedToZero = overloaded{
    core::element::arg_list<std::tuple<double, bool>, std::tuple<float, bool>,
                            std::tuple<bool, bool>, std::tuple<int64_t, bool>,
                            std::tuple<int32_t, bool>>,
    [](const auto &a, bool isMasked) { return isMasked ? decltype(a){0} : a; },
    [](const scipp::units::Unit &a, const scipp::units::Unit &b) {
      if (b != scipp::units::dimensionless) {
        throw except::UnitError("Expected mask to contain dimensionless units");
      }

      return a;
    }};

/// Set the elements referenced by a span to 0
template <class T> void zero(const scipp::span<T> &data) {
  for (auto &x : data)
    x = 0.0;
}

/// Set the elements referenced by a span to +Inf
template <class T> void inf(const scipp::span<T> &data) {
auto posinf = std::numeric_limits<T>::infinity();
  for (auto &x : data)
    x = posinf;
}

/// Set the elements referenced by a span to -Inf
template <class T> void ninf(const scipp::span<T> &data) {
auto neginf = - std::numeric_limits<T>::infinity();
  for (auto &x : data)
    x = neginf;
}

/// Set the elements references by the spans for values and variances to 0
template <class T> void zero(const core::ValueAndVariance<span<T>> &data) {
  zero(data.value);
  zero(data.variance);
}

/// Set the elements references by the spans for values and variances to +Inf
template <class T> void inf(const core::ValueAndVariance<span<T>> &data) {
  inf(data.value);
  inf(data.variance);
}

/// Set the elements references by the spans for values and variances to -Inf
template <class T> void ninf(const core::ValueAndVariance<span<T>> &data) {
  ninf(data.value);
  ninf(data.variance);
}

} // namespace scipp::core::element
