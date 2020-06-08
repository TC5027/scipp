// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2020 Scipp contributors (https://github.com/scipp)
/// @file
/// @author Simon Heybrock
#pragma once

#include "scipp/dataset/dataset.h"
#include "scipp/variable/rebin.h"

namespace scipp::dataset {

[[nodiscard]] SCIPP_DATASET_EXPORT DataArray rebin(
    const DataArrayConstView &a, const Dim dim, const VariableConstView &coord);
[[nodiscard]] SCIPP_DATASET_EXPORT Dataset
rebin(const DatasetConstView &d, const Dim dim, const VariableConstView &coord);

[[nodiscard]] SCIPP_DATASET_EXPORT DataArray resample(const DataArrayConstView &a, const Dim dim,
                                     const VariableConstView &coord,
                                     const scipp::variable::ResampleOp op = scipp::variable::ResampleOp::Sum);
[[nodiscard]] SCIPP_DATASET_EXPORT Dataset resample(const DatasetConstView &d, const Dim dim,
                                   const VariableConstView &coord,
                                   const scipp::variable::ResampleOp op = scipp::variable::ResampleOp::Sum);

} // namespace scipp::dataset
