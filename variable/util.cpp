// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2020 Scipp contributors (https://github.com/scipp)
/// @file
/// @author Simon Heybrock
#include "scipp/variable/util.h"
#include "scipp/core/element/util.h"
#include "scipp/core/except.h"
#include "scipp/variable/arithmetic.h"
#include "scipp/variable/comparison.h"
#include "scipp/variable/except.h"
#include "scipp/variable/misc_operations.h"
#include "scipp/variable/reduction.h"
#include "scipp/variable/transform.h"

using namespace scipp::core;

namespace scipp::variable {

Variable linspace(const VariableConstView &start, const VariableConstView &stop,
                  const Dim dim, const scipp::index num) {
  // The implementation here is slightly verbose and explicit. It could be
  // improved if we were to introduce new variants of `transform`, similar to
  // `std::generate`.
  core::expect::equals(start.dims(), stop.dims());
  core::expect::equals(start.unit(), stop.unit());
  core::expect::equals(start.dtype(), stop.dtype());
  if (start.dtype() != dtype<double> && start.dtype() != dtype<float>)
    throw except::TypeError(
        "Cannot create linspace with non-floating-point start and/or stop.");
  if (start.hasVariances() || stop.hasVariances())
    throw except::VariancesError(
        "Cannot create linspace with start and/or stop containing variances.");
  auto dims = start.dims();
  dims.addInner(dim, num);
  Variable out(start, dims);
  const auto range = stop - start;
  for (scipp::index i = 0; i < num - 1; ++i)
    out.slice({dim, i}).assign(
        start +
        astype(static_cast<double>(i) / (num - 1) * units::one, start.dtype()) *
            range);
  out.slice({dim, num - 1}).assign(stop); // endpoint included
  return out;
}

Variable values(const VariableConstView &x) {
  return transform(x, element::values);
}
Variable variances(const VariableConstView &x) {
  return transform(x, element::variances);
}

/// Return true if variable values are sorted along given dim.
///
/// If `order` is SortOrder::Ascending, checks if values are non-decreasing.
/// If `order` is SortOrder::Descending, checks if values are non-increasing.
bool is_sorted(const VariableConstView &x, const Dim dim,
               const SortOrder order) {
  const auto size = x.dims()[dim];
  if (size < 2)
    return true;
  auto out = makeVariable<bool>(Values{true});
  if (order == SortOrder::Ascending)
    accumulate_in_place(out, x.slice({dim, 0, size - 1}),
                        x.slice({dim, 1, size}),
                        core::element::is_sorted_nondescending);
  else
    accumulate_in_place(out, x.slice({dim, 0, size - 1}),
                        x.slice({dim, 1, size}),
                        core::element::is_sorted_nonascending);
  return out.value<bool>();
}

} // namespace scipp::variable
