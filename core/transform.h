// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2019 Scipp contributors (https://github.com/scipp)
/// @file
/// @author Simon Heybrock
#ifndef TRANSFORM_H
#define TRANSFORM_H

#include "except.h"
#include "variable.h"
#include "visit.h"

namespace scipp::core {

namespace detail {

template <class T> struct ValueAndVariance {
  T value;
  T variance;
};

template <class T1, class T2>
constexpr auto operator+(const ValueAndVariance<T1> a,
                         const ValueAndVariance<T2> b) noexcept {
  return ValueAndVariance{a.value + b.value, a.variance + b.variance};
}
template <class T1, class T2>
constexpr auto operator-(const ValueAndVariance<T1> a,
                         const ValueAndVariance<T2> b) noexcept {
  return ValueAndVariance{a.value - b.value, a.variance - b.variance};
}
template <class T1, class T2>
constexpr auto operator*(const ValueAndVariance<T1> a,
                         const ValueAndVariance<T2> b) noexcept {
  return ValueAndVariance{a.value * b.value,
                          a.variance * b.value * b.value +
                              b.variance * a.value * a.value};
}
template <class T1, class T2>
constexpr auto operator/(const ValueAndVariance<T1> a,
                         const ValueAndVariance<T2> b) noexcept {
  return ValueAndVariance{
      a.value / b.value,
      (a.variance + b.variance * (a.value * a.value) / (b.value * b.value)) /
          (b.value * b.value)};
}

template <class T1, class T2>
constexpr auto operator+(const ValueAndVariance<T1> a, const T2 b) noexcept {
  return ValueAndVariance{a.value + b, a.variance};
}
template <class T1, class T2>
constexpr auto operator-(const ValueAndVariance<T1> a, const T2 b) noexcept {
  return ValueAndVariance{a.value - b, a.variance};
}
template <class T1, class T2>
constexpr auto operator*(const ValueAndVariance<T1> a, const T2 b) noexcept {
  return ValueAndVariance{a.value * b, a.variance * b * b};
}
template <class T1, class T2>
constexpr auto operator/(const ValueAndVariance<T1> a, const T2 b) noexcept {
  return ValueAndVariance{a.value / b, a.variance / (b * b)};
}

template <class T>
ValueAndVariance(const T &val, const T &var)->ValueAndVariance<T>;

template <class T>
std::unique_ptr<VariableConceptT<T>>
makeVariableConceptT(const Dimensions &dims);
template <class T>
std::unique_ptr<VariableConceptT<T>>
makeVariableConceptT(const Dimensions &dims, Vector<T> data);

template <class Op> struct TransformSparse {
  Op op;
  // TODO avoid copies... need in place transform (for_each, but with a second
  // input range).
  template <class T> constexpr auto operator()(sparse_container<T> x) const {
    std::transform(x.begin(), x.end(), x.begin(), op);
    return x;
  }
  // TODO Would like to use T1 and T2 for a and b, but currently this leads to
  // selection of the wrong overloads.
  template <class T>
  constexpr auto operator()(sparse_container<T> a, const T b) const {
    std::transform(a.begin(), a.end(), a.begin(),
                   [&, b](const T a) { return op(a, b); });
    return a;
  }
};

template <class T> struct is_eigen_type : std::false_type {};
template <class T, int Rows, int Cols>
struct is_eigen_type<Eigen::Matrix<T, Rows, Cols>> : std::true_type {};
template <class T>
inline constexpr bool is_eigen_type_v = is_eigen_type<T>::value;

template <class T> struct is_sparse : std::false_type {};
template <class T> struct is_sparse<sparse_container<T>> : std::true_type {};
template <class T> inline constexpr bool is_sparse_v = is_sparse<T>::value;

template <class T> struct as_view {
  using value_type = typename T::value_type;
  bool hasVariances() const { return data.hasVariances(); }
  auto values() const { return data.valuesView(dims); }
  auto variances() const { return data.variancesView(dims); }

  T &data;
  const Dimensions &dims;
};
template <class T> as_view(T &data, const Dimensions &dims)->as_view<T>;

template <class T1, class T2, class Op>
void do_transform(const T1 &a, const T2 &b, T1 &c, Op op) {
  auto a_val = a.values();
  auto b_val = b.values();
  auto c_val = c.values();
  if (a.hasVariances()) {
    if constexpr (is_sparse_v<typename T1::value_type> ||
                  is_sparse_v<typename T2::value_type>) {
      throw std::runtime_error(
          "Propagation of uncertainties for sparse data not implemented yet.");
    } else if constexpr (is_eigen_type_v<typename T1::value_type> ||
                         is_eigen_type_v<typename T2::value_type>) {
      throw std::runtime_error("This dtype cannot have a variance.");
    } else {
      auto a_var = a.variances();
      auto c_var = c.variances();
      if (b.hasVariances()) {
        auto b_var = b.variances();
        for (scipp::index i = 0; i < a_val.size(); ++i) {
          const ValueAndVariance a_{a_val[i], a_var[i]};
          const ValueAndVariance b_{b_val[i], b_var[i]};
          const auto out = op(a_, b_);
          c_val[i] = out.value;
          c_var[i] = out.variance;
        }
      } else {
        for (scipp::index i = 0; i < a_val.size(); ++i) {
          const ValueAndVariance a_{a_val[i], a_var[i]};
          const auto out = op(a_, b_val[i]);
          c_val[i] = out.value;
          c_var[i] = out.variance;
        }
      }
    }
  } else if (b.hasVariances()) {
    throw std::runtime_error(
        "RHS in operation has variances but LHS does not.");
  } else {
    std::transform(a_val.begin(), a_val.end(), b_val.begin(), c_val.begin(),
                   op);
  }
}

template <class Op> struct TransformInPlace {
  Op op;
  template <class T> void operator()(T &&handle) const {
    auto data = handle->values();
    std::transform(data.begin(), data.end(), data.begin(), op);
  }
  template <class A, class B> void operator()(A &&a, B &&b_ptr) const {
    // std::unique_ptr::operator*() is const but returns mutable reference, need
    // to artificially put const to we call the correct overloads of ViewModel.
    const auto &b = *b_ptr;
    const auto &dimsA = a->dims();
    const auto &dimsB = b.dims();
    try {
      if constexpr (std::is_same_v<decltype(*a), decltype(*b_ptr)>) {
        if (a->valuesView(dimsA).overlaps(b.valuesView(dimsA))) {
          // If there is an overlap between lhs and rhs we copy the rhs before
          // applying the operation.
          const auto &data = b.valuesView(b.dims());
          using T = typename std::remove_reference_t<decltype(b)>::value_type;
          const std::unique_ptr<VariableConceptT<T>> copy =
              detail::makeVariableConceptT<T>(
                  dimsB, Vector<T>(data.begin(), data.end()));
          return operator()(a, copy);
        }
      }

      if (a->isContiguous() && dimsA.contains(dimsB)) {
        if (b.isContiguous() && dimsA.isContiguousIn(dimsB)) {
          do_transform(*a, b, *a, op);
        } else {
          do_transform(*a, as_view{b, dimsA}, *a, op);
        }
      } else if (dimsA.contains(dimsB)) {
        auto a_view = as_view{*a, dimsA};
        if (b.isContiguous() && dimsA.isContiguousIn(dimsB)) {
          do_transform(a_view, b, a_view, op);
        } else {
          do_transform(a_view, as_view{b, dimsA}, a_view, op);
        }
      } else {
        // LHS has fewer dimensions than RHS, e.g., for computing sum. Use view.
        auto a_view = as_view{*a, dimsB};
        if (b.isContiguous() && dimsA.isContiguousIn(dimsB)) {
          do_transform(a_view, b, a_view, op);
        } else {
          do_transform(a_view, as_view{b, dimsB}, a_view, op);
        }
      }
    } catch (const std::bad_cast &) {
      throw std::runtime_error("Cannot apply arithmetic operation to "
                               "Variables: Underlying data types do not "
                               "match.");
    }
  }
};
template <class Op> TransformInPlace(Op)->TransformInPlace<Op>;

template <class Op> struct Transform {
  Op op;
  template <class T> VariableConceptHandle operator()(T &&handle) const {
    if (handle->hasVariances())
      throw std::runtime_error(
          "Propgation of uncertainties not implemented for this case.");
    auto data = handle->values();
    // TODO Should just make empty container here, without init.
    auto out = detail::makeVariableConceptT<decltype(op(*data.begin()))>(
        handle->dims());
    // TODO Typo data->values() also compiles, but const-correctness should
    // forbid this.
    auto outData = out->values();
    std::transform(data.begin(), data.end(), outData.begin(), op);
    return {std::move(out)};
  }
};

} // namespace detail

template <class... Ts> struct overloaded : Ts... { using Ts::operator()...; };
template <class... Ts> overloaded(Ts...)->overloaded<Ts...>;

/// Transform the data elements of a variable in-place.
//
// Note that this is deliberately not named `for_each`: Unlike std::for_each,
// this function does not promise in-order execution. This overload is
// equivalent to std::transform with a single input range and an output range
// identical to the input range, but avoids potentially costly element copies.
template <class... Ts, class Var, class Op>
void transform_in_place(Var &var, Op op) {
  using namespace detail;
  try {
    scipp::core::visit_impl<Ts...>::apply(
        TransformInPlace{overloaded{op, TransformSparse<Op>{op}}},
        var.dataHandle().variant());
  } catch (const std::bad_variant_access &) {
    throw std::runtime_error("Operation not implemented for this type.");
  }
}

/// Transform the data elements of a variable in-place.
//
// This overload is equivalent to std::transform with two input ranges and an
// output range identical to the secound input range, but avoids potentially
// costly element copies.
template <class... TypePairs, class Var1, class Var, class Op>
void transform_in_place(const Var1 &other, Var &&var, Op op) {
  using namespace detail;
  try {
    scipp::core::visit(std::tuple_cat(TypePairs{}...))
        .apply(TransformInPlace{overloaded{op, TransformSparse<Op>{op}}},
               var.dataHandle().variant(), other.dataHandle().variant());
  } catch (const std::bad_variant_access &) {
    throw except::TypeError("Cannot apply operation to item dtypes " +
                            to_string(var.dtype()) + " and " +
                            to_string(other.dtype()) + '.');
  }
}

/// Transform the data elements of a variable and return a new Variable.
//
// This overload is equivalent to std::transform with a single input range, but
// avoids the need to manually create a new variable for the output and the need
// for, e.g., std::back_inserter.
template <class... Ts, class Var, class Op>
Variable transform(const Var &var, Op op) {
  using namespace detail;
  try {
    return Variable(var, scipp::core::visit_impl<Ts...>::apply(
                             Transform<Op>{op}, var.dataHandle().variant()));
  } catch (const std::bad_variant_access &) {
    throw std::runtime_error("Operation not implemented for this type.");
  }
}

} // namespace scipp::core

#endif // TRANSFORM_H