// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2020 Scipp contributors (https://github.com/scipp)
/// @file
/// @author Simon Heybrock
#pragma once

#include "scipp/core/dimensions.h"

namespace scipp::variable {

bool isBinEdge(const Dim dim, Dimensions edges, const Dimensions &toMatch);

} // namespace scipp::variable
