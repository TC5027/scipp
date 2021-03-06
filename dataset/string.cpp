// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2020 Scipp contributors (https://github.com/scipp)
/// @file
/// @author Simon Heybrock
#include <algorithm> // for std::sort
#include <set>

#include "scipp/dataset/dataset.h"
#include "scipp/dataset/except.h"

namespace scipp::dataset {

std::ostream &operator<<(std::ostream &os, const DataArrayConstView &data) {
  return os << to_string(data);
}

std::ostream &operator<<(std::ostream &os, const DataArrayView &data) {
  return os << DataArrayConstView(data);
}

std::ostream &operator<<(std::ostream &os, const DataArray &data) {
  return os << DataArrayConstView(data);
}

std::ostream &operator<<(std::ostream &os, const DatasetConstView &dataset) {
  return os << to_string(dataset);
}

std::ostream &operator<<(std::ostream &os, const DatasetView &dataset) {
  return os << DatasetConstView(dataset);
}

std::ostream &operator<<(std::ostream &os, const Dataset &dataset) {
  return os << DatasetConstView(dataset);
}

constexpr const char *tab = "    ";

template <class D>
std::string do_to_string(const D &dataset, const std::string &id,
                         const Dimensions &dims, const std::string &shift = "");

template <class T> auto sorted(const T &map) {
  using core::to_string;
  std::vector<std::pair<std::string, VariableConstView>> elems;
  for (const auto &[dim, var] : map)
    elems.emplace_back(to_string(dim), var);
  std::sort(elems.begin(), elems.end(),
            [](const auto &a, const auto &b) { return a.first < b.first; });
  return elems;
}

template <class Key>
auto format_data_view(const Key &name, const DataArrayConstView &data,
                      const Dimensions &datasetDims = Dimensions()) {
  std::stringstream s;
  if (data.hasData())
    s << format_variable(name, data.data(), datasetDims);
  else {
    s << tab << name << " (data not histogrammed yet)\n";
    s << tab << "Unaligned:\n";
    s << do_to_string(data.unaligned(), "", data.unaligned().dims(),
                      std::string(tab) + tab);
  }
  if (!data.masks().empty()) {
    s << tab << "Masks:\n";
    for (const auto &[key, var] : sorted(data.masks()))
      s << tab << tab << format_variable(key, var, datasetDims);
  }
  if (!data.unaligned_coords().empty()) {
    s << tab << "Coordinates (unaligned):\n";
    for (const auto &[key, var] : sorted(data.unaligned_coords()))
      s << tab << tab << format_variable(key, var, datasetDims);
  }
  return s.str();
}

template <class D>
std::string do_to_string(const D &dataset, const std::string &id,
                         const Dimensions &dims, const std::string &shift) {
  std::stringstream s;
  if (!id.empty())
    s << shift << id + '\n';
  s << shift << "Dimensions: " << to_string(dims) << '\n';

  if (!dataset.coords().empty()) {
    s << shift << "Coordinates:\n";
    CoordsConstView map;
    if constexpr (std::is_same_v<D, DataArray> ||
                  std::is_same_v<D, DataArrayConstView>)
      map = dataset.aligned_coords();
    else
      map = dataset.coords();
    for (const auto &[name, var] : sorted(map))
      s << shift << format_variable(name, var, dims);
  }

  if constexpr (std::is_same_v<D, DataArray> ||
                std::is_same_v<D, DataArrayConstView>) {
    s << shift << "Data:\n" << format_data_view(dataset.name(), dataset);
  } else {
    if (!dataset.empty())
      s << shift << "Data:\n";
    std::set<std::string> sorted_items;
    for (const auto &item : dataset)
      sorted_items.insert(item.name());
    for (const auto &name : sorted_items)
      s << shift << format_data_view(name, dataset[name], dims);
  }

  s << '\n';
  return s.str();
}

template <class T> Dimensions dimensions(const T &dataset) {
  Dimensions dims;
  for (const auto &[dim, size] : dataset.dimensions())
    dims.add(dim, size);
  return dims;
}

std::string to_string(const DataArray &data) {
  return do_to_string(data, "<scipp.DataArray>", data.dims());
}

std::string to_string(const DataArrayConstView &data) {
  return do_to_string(data, "<scipp.DataArrayView>", data.dims());
}

std::string to_string(const Dataset &dataset) {
  return do_to_string(dataset, "<scipp.Dataset>", dimensions(dataset));
}

std::string to_string(const DatasetConstView &dataset) {
  return do_to_string(dataset, "<scipp.DatasetView>", dimensions(dataset));
}

} // namespace scipp::dataset
