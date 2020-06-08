// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2020 Scipp contributors (https://github.com/scipp)
/// @file
/// @author Simon Heybrock
#pragma once

#include "scipp/core/dimensions.h"
#include "scipp-variable_export.h"
#include "scipp/variable/variable.h"

namespace scipp::variable {

[[nodiscard]] SCIPP_VARIABLE_EXPORT Variable
rebin(const VariableConstView &var, const Dim dim,
      const VariableConstView &oldCoord, const VariableConstView &newCoord);

enum class ResampleOp { Sum, Mean, Min, Max };
[[nodiscard]] SCIPP_VARIABLE_EXPORT Variable resample(const VariableConstView &var,
                                     const Dim dim,
                                     const VariableConstView &oldCoord,
                                     const VariableConstView &newCoord,
                                     const ResampleOp op = ResampleOp::Sum);


bool isBinEdge(const Dim dim, Dimensions edges, const Dimensions &toMatch);

} // namespace scipp::variable
