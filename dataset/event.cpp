// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2020 Scipp contributors (https://github.com/scipp)
/// @file
/// @author Simon Heybrock
#include <limits>

#include "scipp/core/event.h"
#include "scipp/core/subspan_view.h"
#include "scipp/core/transform.h"
#include "scipp/core/variable_operations.h"

#include "scipp/dataset/dataset.h"
#include "scipp/dataset/event.h"

namespace scipp::dataset {
/// Return true if a data array contains events
bool is_events(const DataArrayConstView &array) {
  if (array.hasData() && is_events(array.data()))
    return true;
  for (const auto &item : array.coords())
    if (is_events(item.second))
      return true;
  return false;
}

namespace event {

void append(const DataArrayView &a, const DataArrayConstView &b) {
  if (!is_events(a) || !is_events(b))
    throw except::EventDataError("Cannot concatenate non-event data.");

  if (is_events(a.data())) {
    core::event::append(a.data(),
                        is_events(b.data()) ? b.data() : broadcast_weights(b));
  } else if (is_events(b.data())) {
    a.setData(core::event::concatenate(broadcast_weights(a), b.data()));
  } else if (a.data() != b.data()) {
    a.setData(
        core::event::concatenate(broadcast_weights(a), broadcast_weights(b)));
  } else {
    // Do nothing for identical scalar weights
  }
  for (const auto &[dim, coord] : a.coords())
    if (is_events(coord))
      core::event::append(coord, b.coords()[dim]);
    else
      core::expect::equals(coord, b.coords()[dim]);
}

DataArray concatenate(const DataArrayConstView &a,
                      const DataArrayConstView &b) {
  DataArray out(a);
  append(out, b);
  return out;
}

/// Broadcast scalar weights of data array containing event data.
Variable broadcast_weights(const DataArrayConstView &events) {
  for (const auto &item : events.coords())
    if (is_events(item.second))
      return core::event::broadcast(events.data(), item.second);
  throw except::EventDataError(
      "No coord with event lists found, cannot broadcast weights.");
}

namespace filter_detail {
template <class T>
using make_select_args = std::tuple<event_list<T>, span<const T>>;
template <class T, class Index>
using copy_if_args = std::tuple<event_list<T>, event_list<Index>>;

/// Return new variable with values copied from `var` if index is included in
/// `select`.
constexpr auto copy_if = [](const VariableConstView &var,
                            const VariableConstView &select) {
  return core::transform<std::tuple<
      copy_if_args<double, int32_t>, copy_if_args<float, int32_t>,
      copy_if_args<int64_t, int32_t>, copy_if_args<int32_t, int32_t>,
      copy_if_args<double, int64_t>, copy_if_args<float, int64_t>,
      copy_if_args<int64_t, int64_t>, copy_if_args<int32_t, int64_t>>>(
      var, select,
      overloaded{
          core::transform_flags::expect_no_variance_arg<1>,
          [](const auto &var_, const auto &select_) {
            using VarT = std::decay_t<decltype(var_)>;
            using Events = event_list<typename VarT::value_type>;
            const auto size = scipp::size(select_);
            if constexpr (core::detail::is_ValuesAndVariances_v<VarT>) {
              std::pair<Events, Events> out;
              out.first.reserve(size);
              out.second.reserve(size);
              for (const auto i : select_) {
                out.first.push_back(var_.values[i]);
                out.second.push_back(var_.variances[i]);
              }
              return out;
            } else {
              Events out;
              out.reserve(size);
              for (const auto i : select_)
                out.push_back(var_[i]);
              return out;
            }
          },
          [](const units::Unit &var_, const units::Unit &) { return var_; }});
};

/// Return list of indices with coord values for given dim inside interval.
template <class T>
const auto make_select = [](const DataArrayConstView &array, const Dim dim,
                            const VariableConstView &interval) {
  return core::transform<
      std::tuple<make_select_args<double>, make_select_args<float>,
                 make_select_args<int64_t>, make_select_args<int32_t>>>(
      array.coords()[dim], subspan_view(interval, dim),
      overloaded{core::transform_flags::expect_no_variance_arg<0>,
                 core::transform_flags::expect_no_variance_arg<1>,
                 [](const auto &coord_, const auto &interval_) {
                   const auto low = interval_[0];
                   const auto high = interval_[1];
                   const auto size = scipp::size(coord_);
                   event_list<T> select_;
                   for (scipp::index i = 0; i < size; ++i)
                     if (coord_[i] >= low && coord_[i] < high)
                       select_.push_back(i);
                   return select_;
                 },
                 [](const units::Unit &coord_, const units::Unit &interval_) {
                   core::expect::equals(coord_, interval_);
                   return units::Unit(units::dimensionless);
                 }});
};

} // namespace filter_detail

/// Return filtered event data based on excluding all events with coord values
/// for given dim outside interval.
DataArray filter(const DataArrayConstView &array, const Dim dim,
                 const VariableConstView &interval,
                 const AttrPolicy attrPolicy) {
  using namespace filter_detail;
  const auto &max_event_list_length =
      max(core::event::sizes(array.coords()[dim]));
  const bool need_64bit_indices =
      max_event_list_length.values<scipp::index>()[0] >
      std::numeric_limits<int32_t>::max();
  const auto select = need_64bit_indices
                          ? make_select<int64_t>(array, dim, interval)
                          : make_select<int32_t>(array, dim, interval);

  std::map<Dim, Variable> coords;
  for (const auto &[d, coord] : array.coords())
    coords.emplace(d, is_events(coord) ? copy_if(coord, select) : copy(coord));

  Dataset empty;
  return DataArray{is_events(array.data()) ? copy_if(array.data(), select)
                                           : copy(array.data()),
                   std::move(coords), array.masks(),
                   attrPolicy == AttrPolicy::Keep ? array.attrs()
                                                  : empty.attrs()};
}

} // namespace event
} // namespace scipp::dataset