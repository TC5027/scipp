// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2020 Scipp contributors (https://github.com/scipp)
/// @file
/// @author Simon Heybrock
#include "scipp/dataset/rebin.h"
#include "scipp/common/overloaded.h"
#include "scipp/variable/rebin.h"

#include "dataset_operations_common.h"

namespace scipp::dataset {

DataArray rebin(const DataArrayConstView &a, const Dim dim,
                const VariableConstView &coord) {
  auto rebinned = apply_to_data_and_drop_dim(
      a,
      overloaded{no_realigned_support,
                 [](auto &&... _) { return rebin(_...); }},
      dim, a.coords()[dim], coord);

  for (auto &&[name, mask] : a.masks()) {
    if (mask.dims().contains(dim))
      rebinned.masks().set(name, rebin(mask, dim, a.coords()[dim], coord));
  }

  rebinned.coords().set(dim, coord);
  return rebinned;
}

Dataset rebin(const DatasetConstView &d, const Dim dim,
              const VariableConstView &coord) {
  return apply_to_items(
      d, [](auto &&... _) { return rebin(_...); }, dim, coord);
}

DataArray resample(const DataArrayConstView &a, const Dim dim,
                const VariableConstView &coord, const scipp::variable::ResampleOp op) {
  auto resampled = apply_to_data_and_drop_dim(
      a,
      overloaded{no_realigned_support,
                 [](auto &&... _) { return resample(_...); }},
      dim, a.coords()[dim], coord, op);

  for (auto &&[name, mask] : a.masks()) {
    if (mask.dims().contains(dim))
      resampled.masks().set(name, resample(mask, dim, a.coords()[dim], coord, op));
  }

  resampled.coords().set(dim, coord);
  return resampled;
}

Dataset resample(const DatasetConstView &d, const Dim dim,
              const VariableConstView &coord, const scipp::variable::ResampleOp op) {
  return apply_to_items(
      d, [](auto &&... _) { return resample(_...); }, dim, coord, op);
}

} // namespace scipp::dataset
