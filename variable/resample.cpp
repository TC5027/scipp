// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2020 Scipp contributors (https://github.com/scipp)
/// @file
/// @author Simon Heybrock, Igor Gudich
#include "scipp/core/element/resample.h"
#include "scipp/units/except.h"
#include "scipp/variable/apply.h"
#include "scipp/variable/arithmetic.h"
#include "scipp/variable/except.h"
#include "scipp/variable/misc_operations.h"
#include "scipp/variable/reduction.h"
#include "scipp/variable/rebin.h"
#include "scipp/variable/transform_subspan.h"
#include <iostream>

namespace scipp::variable {

// bool isBinEdge(const Dim dim, Dimensions edges, const Dimensions &toMatch) {
//   edges.resize(dim, edges[dim] - 1);
//   return edges[dim] == toMatch[dim];
// }

template <typename T, typename Lambda>
void do_resample_non_inner(const Dim dim, const VariableConstView &oldT,
                     Variable &newT, const VariableConstView &oldCoordT,
                     const VariableConstView &newCoordT, Lambda &&lambda) {
  const auto oldSize = oldT.dims()[dim];
  const auto newSize = newT.dims()[dim];
  // Variable counter{makeVariable<T>(Dims{dim}, Shape{newSize})};
  // std::vector<int> newcounter(newT.dims()[dim]);

  const auto *xold = oldCoordT.values<T>().data();
  const auto *xnew = newCoordT.values<T>().data();
  // This function assumes that dimensions between coord and data
  // coord is 1D.
  int iold = 0;
  int inew = 0;
  while ((iold < oldSize) && (inew < newSize)) {
    auto xo_low = xold[iold];
    auto xo_high = xold[iold + 1];
    auto xn_low = xnew[inew];
    auto xn_high = xnew[inew + 1];

    if (xn_high <= xo_low)
      inew++; /* old and new bins do not overlap */
    else if (xo_high <= xn_low)
      iold++; /* old and new bins do not overlap */
    else {
      // delta is the overlap of the bins on the x axis
      // auto delta = std::min(xn_high, xo_high) - std::max(xn_low, xo_low);

      // auto owidth = xo_high - xo_low;
 

      // auto newvals = max(concatenate(newT.slice({dim, inew}),
      //     astype(oldT.slice({dim, iold}),
      //            newT.dtype()), dim), dim);
      // newT.slice({dim, inew}) -= newT.slice({dim, inew});
      // newT.slice({dim, inew}) += newvals;
 
      // newT.slice({dim, inew}) = newT.slice({dim, inew}) / counter.slice({dim, inew}) 

      // // Implement running mean
      // if (counter.slice({dim, inew}).values<T>()[0] > 0.0)
      //    newT.slice({dim, inew}) /= counter.slice({dim, inew});
      // newT.slice({dim, inew}) *= counter.slice({dim, inew}) / (counter.slice({dim, inew}) + 1.0 * units::one);
      // counter.slice({dim, inew}) += 1.0 * units::one;
      // newT.slice({dim, inew}) += astype(oldT.slice({dim, iold}),
      //            newT.dtype()) / counter.slice({dim, inew});

      // Sum implementation
      // auto newvals = max(concatenate(newT.slice({dim, inew}),
      //     astype(oldT.slice({dim, iold}),
      //            newT.dtype()), dim), dim);
      // newT.slice({dim, inew}) -= newT.slice({dim, inew});

      // lambda(newT.slice({dim, inew}),astype(oldT.slice({dim, iold}),
      //            newT.dtype()));
      lambda(newT, oldT, dim, inew, iold);


      // counter.slice({dim, inew}) += 1.0 * units::one;


      // newT.slice({dim, inew}) = max(concatenate(newT.slice({dim, inew}),
      //     // astype(oldT.slice({dim, iold}) * ((delta / owidth) * units::one),
      //     astype(oldT.slice({dim, iold}),
      //            newT.dtype()), dim), dim);
      if (xn_high > xo_high) {
        iold++;
      } else {
        inew++;
      }
    }
  }
}

// static constexpr auto resample_sum = [](Variable &newT, const VariableConstView &oldT, const Dim dim, const int inew, const int iold){ newT.slice({dim, inew}) += astype(oldT.slice({dim, iold}),
//                  newT.dtype());};

template <typename T>
void resample_non_inner_sum(const Dim dim, const VariableConstView &oldT,
                     Variable &newT, const VariableConstView &oldCoordT,
                     const VariableConstView &newCoordT) {
  constexpr auto lambda = [](Variable &newT_, const VariableConstView &oldT_, const Dim dim_, const int inew, const int iold){ newT_.slice({dim_, inew}) += astype(oldT_.slice({dim_, iold}),
                 newT_.dtype());};
  return do_resample_non_inner<T>(dim, oldT, newT, oldCoordT, newCoordT, lambda);
}


// static constexpr auto resample_max = [](Variable &newT, const VariableConstView &oldT, const Dim dim, const int inew, const int iold){ newT.slice({dim, inew}) += astype(oldT.slice({dim, iold}),
//                  newT.dtype());};

template <typename T>
void resample_non_inner_max(const Dim dim, const VariableConstView &oldT,
                     Variable &newT, const VariableConstView &oldCoordT,
                     const VariableConstView &newCoordT) {

  constexpr auto lambda = [](Variable &newT_, const VariableConstView &oldT_, const Dim dim_, const int inew, const int iold){
      auto newvals = max(concatenate(newT_.slice({dim_, inew}),
          astype(oldT_.slice({dim_, iold}),
                 newT_.dtype()), dim_), dim_);
      // newT.slice({dim, inew}) -= newT.slice({dim, inew});
      newT_.slice({dim_, inew}) += newvals - newT_.slice({dim_, inew});};

    // newT.slice({dim, inew}) += astype(oldT.slice({dim, iold}),
    //              newT.dtype());};

  return do_resample_non_inner<T>(dim, oldT, newT, oldCoordT, newCoordT, lambda);
}


template <class T>
void resample_non_inner(const Dim dim, const VariableConstView &var, 
                Variable &resampled,
               const VariableConstView &oldCoord,
               const VariableConstView &newCoord,
               const ResampleOp op) {
      if (op == ResampleOp::Sum)
        resample_non_inner_sum<double>(dim, var, resampled, oldCoord, newCoord);
      else if (op == ResampleOp::Max)
        resample_non_inner_max<double>(dim, var, resampled, oldCoord, newCoord);
      else
        throw std::runtime_error(
          "Unknown resammpling operation.");
}

namespace resample_inner_detail {
template <class Out, class OutEdge, class In, class InEdge>
using args = std::tuple<span<Out>, span<const OutEdge>, span<const In>,
                        span<const InEdge>>;
}

Variable resample(const VariableConstView &var, const Dim dim,
               const VariableConstView &oldCoord,
               const VariableConstView &newCoord,
               const ResampleOp op) {
  // Rebin could also implemented for count-densities. However, it may be better
  // to avoid this since it increases complexity. Instead, densities could
  // always be computed on-the-fly for visualization, if required.
  // core::expect::unit_any_of(var, {units::counts, units::one});
  if (!isBinEdge(dim, oldCoord.dims(), var.dims()))
    throw except::BinEdgeError(
        "The input does not have coordinates with bin-edges.");

  if (var.dims().inner() == dim) {
    // std::cout << "Using transform_subspan" << std::endl;
    using namespace resample_inner_detail;
    // const auto newSize = newT.dims()[dim];
    // Variable counter{newCoord.slice({dim, 0, newCoord.dims()[dim] - 1})};
    // counter.setUnit(units::one);
    // zero(counter);
    return transform_subspan<std::tuple<
        args<double, double, double, double>, args<float, float, float, float>,
        args<float, double, float, double>, args<float, float, float, double>,
        args<bool, double, bool, double>>>(var.dtype(), dim,
                                           newCoord.dims()[dim] - 1, newCoord,
                                           var, oldCoord, core::element::resample);

  } else {
    // std::cout << "Using rebin_non_inner" << std::endl;
    auto dims = var.dims();
    dims.resize(dim, newCoord.dims()[dim] - 1);
    Variable resampled(var, dims);
    // Variable counts = var / var;
    if (newCoord.dims().ndim() > 1)
      throw std::runtime_error(
          "Not inner rebin works only for 1d coordinates for now.");
    if (oldCoord.dtype() == dtype<double>)
      resample_non_inner<double>(dim, var, resampled, oldCoord, newCoord, op);
      // if (op == ResampleOp::Sum)
      //   resample_non_inner_sum<double>(dim, var, resampled, oldCoord, newCoord);
      // else if (op == ResampleOp::Max)
      //   resample_non_inner_max<double>(dim, var, resampled, oldCoord, newCoord);
      // else
      //   throw std::runtime_error(
      //     "Unknown resammpling operation.");
    else if (oldCoord.dtype() == dtype<float>)
      resample_non_inner<float>(dim, var, resampled, oldCoord, newCoord, op);
      // if (op == ResampleOp::Sum)
      //   resample_non_inner_sum<float>(dim, var, resampled, oldCoord, newCoord);
      // else if (op == ResampleOp::Max)
      //   resample_non_inner_max<float>(dim, var, resampled, oldCoord, newCoord);
    else
      throw std::runtime_error(
          "Rebinning is possible only for double and float types.");
    return resampled;
  }
}

} // namespace scipp::variable
