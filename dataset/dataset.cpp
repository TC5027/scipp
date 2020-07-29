// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2020 Scipp contributors (https://github.com/scipp)
/// @file
/// @author Simon Heybrock
#include "scipp/dataset/dataset.h"
#include "scipp/common/index.h"
#include "scipp/core/except.h"
#include "scipp/dataset/except.h"
#include "scipp/dataset/unaligned.h"

#include "dataset_operations_common.h"

namespace scipp::dataset {

template <class T> typename T::view_type makeViewItem(T &variable) {
  if constexpr (std::is_const_v<T>)
    return typename T::view_type(typename T::const_view_type(variable));
  else
    return typename T::view_type(variable);
}

namespace {
constexpr auto contains = [](const auto &dims, const Dim dim) {
  if constexpr (std::is_same_v<std::decay_t<decltype(dims)>, Dimensions>)
    return dims.contains(dim);
  else
    return dims.count(dim) == 1;
};
}

template <class T1> auto makeViewItems(T1 &coords) {
  std::unordered_map<typename T1::key_type, VariableView> items;
  for (auto &item : coords)
    items.emplace(item.first, makeViewItem(item.second));
  return items;
}

template <class Dims, class T1>
auto makeViewItems(const Dims &dims, T1 &coords) {
  std::unordered_map<typename T1::key_type, VariableView> items;
  // We preserve only items that are part of the space spanned by the
  // provided parent dimensions.
  auto contained = [&dims](const auto &item2) {
    for (const Dim dim : item2.second.dims().labels())
      if (!contains(dims, dim))
        return false;
    return true;
  };
  for (auto &item : coords)
    if (contained(item))
      items.emplace(item.first, makeViewItem(item.second));
  return items;
}

template <class Item, class Dims, class Attrs>
void addOrthogonalAttrs(Item &items, const Dims &parentDims, Attrs &attrs) {
  for (auto &&[name, attr] : attrs)
    if (!parentDims.contains(attr.dims()))
      items.emplace(name, makeViewItem(attr));
}

template <class Item, class Dims, class Coords>
void addAttrsFromCoords(Item &items, const Dims &parentDims, const Dims &dims,
                        Coords &coords) {
  for (auto &&[dim, coord] : coords) {
    const auto &coordDims = coord.dims();
    if (!coordDims.empty() && contains(parentDims, coordDims.inner()) &&
        !contains(dims, coordDims.inner()))
      items.emplace(dim.name(), makeViewItem(coord));
  }
}

/// Add masks that have become unaligned by slicing.
template <class Item, class Dims, class Coords>
void addMasks(Item &items, const Dims &parentDims, const Dims &dims,
              Coords &coords) {
  for (auto &&[name, coord] : coords) {
    const auto &coordDims = coord.dims();
    if (!coordDims.empty() && contains(parentDims, coordDims.inner()) &&
        !contains(dims, coordDims.inner()))
      items.emplace(name, makeViewItem(coord));
  }
}

Dataset::Dataset(const DatasetConstView &view)
    : Dataset(view, view.coords(), view.masks()) {}

Dataset::Dataset(const DataArrayConstView &data) { setData(data.name(), data); }

Dataset::Dataset(const std::map<std::string, DataArrayConstView> &data) {
  for (const auto &[name, item] : data)
    setData(name, item);
}

/// Removes all data items from the Dataset.
///
/// Coordinates, masks, and attributes are not modified.
/// This operation invalidates any view objects creeated from this dataset.
void Dataset::clear() {
  m_data.clear();
  rebuildDims();
}

/// Return a const view to all coordinates of the dataset.
CoordsConstView Dataset::coords() const noexcept {
  return CoordsConstView(makeViewItems(dimensions(), m_coords));
}

/// Return a view to all coordinates of the dataset.
CoordsView Dataset::coords() noexcept {
  return CoordsView(CoordAccess(this), makeViewItems(dimensions(), m_coords));
}

/// Return a const view to all masks of the dataset.
MasksConstView Dataset::masks() const noexcept {
  return MasksConstView(makeViewItems(dimensions(), m_masks));
}

/// Return a view to all masks of the dataset.
MasksView Dataset::masks() noexcept {
  return MasksView(MaskAccess(this), makeViewItems(dimensions(), m_masks));
}

bool Dataset::contains(const std::string &name) const noexcept {
  return m_data.count(name) == 1;
}

/// Removes a data item from the Dataset
///
/// Coordinates, masks, and attributes are not modified.
/// This operation invalidates any view objects created from this dataset.
void Dataset::erase(const std::string &name) {
  if (m_data.erase(std::string(name)) == 0) {
    throw except::NotFoundError("Expected " + to_string(*this) +
                                " to contain " + name + ".");
  }
  rebuildDims();
}

/// Extract a data item from the Dataset, returning a DataArray
///
/// Aligned coordinates and masks are not modified.
/// This operation invalidates any view objects created from this dataset.
DataArray Dataset::extract(const std::string &name) {
  const auto &view = operator[](name);
  const auto &item = m_data.find(name);

  auto coords = copy_map(view.coords());
  auto masks = copy_map(view.masks());
  auto unaligned_coords = std::move(item->second.coords);

  DataArray extracted;
  if (view.hasData())
    extracted = DataArray(std::move(item->second.data), std::move(coords),
                          std::move(masks), std::move(unaligned_coords), name);
  else
    extracted = DataArray(std::move(*item->second.unaligned), std::move(coords),
                          std::move(masks), std::move(unaligned_coords), name);
  erase(name);
  return extracted;
}

/// Return a const view to data and coordinates with given name.
DataArrayConstView Dataset::operator[](const std::string &name) const {
  scipp::expect::contains(*this, name);
  return DataArrayConstView(*this, *m_data.find(name));
}

/// Return a view to data and coordinates with given name.
DataArrayView Dataset::operator[](const std::string &name) {
  scipp::expect::contains(*this, name);
  return DataArrayView(*this, *m_data.find(name));
}

namespace extents {
// Internally use negative extent -1 to indicate unknown edge state. The `-1`
// is required for dimensions with extent 0.
scipp::index makeUnknownEdgeState(const scipp::index extent) {
  return -extent - 1;
}
scipp::index shrink(const scipp::index extent) { return extent - 1; }
bool isUnknownEdgeState(const scipp::index extent) { return extent < 0; }
scipp::index decodeExtent(const scipp::index extent) {
  if (isUnknownEdgeState(extent))
    return -extent - 1;
  return extent;
}
bool isSame(const scipp::index extent, const scipp::index reference) {
  return reference == -extent - 1;
}
bool oneLarger(const scipp::index extent, const scipp::index reference) {
  return extent == -reference - 1 + 1;
}
bool oneSmaller(const scipp::index extent, const scipp::index reference) {
  return extent == -reference - 1 - 1;
}
void setExtent(std::unordered_map<Dim, scipp::index> &dims, const Dim dim,
               const scipp::index extent, const bool isCoord) {
  const auto it = dims.find(dim);
  // Internally use negative extent -1 to indicate unknown edge state. The `-1`
  // is required for dimensions with extent 0.

  if (it == dims.end()) {
    dims[dim] = extents::makeUnknownEdgeState(extent);
  } else {
    auto &heldExtent = it->second;
    if (extents::isUnknownEdgeState(heldExtent)) {
      if (extents::isSame(extent, heldExtent)) { // Do nothing
      } else if (extents::oneLarger(extent, heldExtent) && isCoord) {
        heldExtent = extents::shrink(extent);
      } else if (extents::oneSmaller(extent, heldExtent) && !isCoord) {
        heldExtent = extent;
      } else {
        throw except::DimensionError(heldExtent, extent);
      }
    } else {
      // Check for known edge state
      if ((extent != heldExtent || isCoord) && extent != heldExtent + 1)
        throw except::DimensionError(heldExtent, extent);
    }
  }
}
} // namespace extents

/// Consistency-enforcing update of the dimensions of the dataset.
///
/// Calling this in the various set* methods prevents insertion of variable with
/// bad shape. This supports insertion of bin edges. Note that the current
/// implementation does not support shape-changing operations which would in
/// theory be permitted but are probably not important in reality: The previous
/// extent of a replaced item is not excluded from the check, so even if that
/// replaced item is the only one in the dataset with that dimension it cannot
/// be "resized" in this way.
void Dataset::setDims(const Dimensions &dims, const Dim coordDim) {
  auto tmp = m_dims;
  for (const auto dim : dims.labels())
    extents::setExtent(tmp, dim, dims[dim], dim == coordDim);
  m_dims = tmp;
}

void Dataset::rebuildDims() {
  m_dims.clear();

  for (const auto &d : *this) {
    setDims(d.dims());
  }
  for (const auto &c : m_coords) {
    setDims(c.second.dims(), dim_of_coord(c.second, c.first));
  }
  for (const auto &m : m_masks) {
    setDims(m.second.dims());
  }
}

/// Set (insert or replace) the coordinate for the given dimension.
void Dataset::setCoord(const Dim dim, Variable coord) {
  setDims(coord.dims(), dim);
  m_coords.insert_or_assign(dim, std::move(coord));
}

/// Set (insert or replace) an unaligned coordinate for item with given name.
void Dataset::setCoord(const std::string &name, const Dim dim, Variable coord) {
  scipp::expect::contains(*this, name);
  setDims(coord.dims(), dim);
  m_data[name].coords.insert_or_assign(dim, std::move(coord));
}

/// Set (insert or replace) the masks for the given mask name.
///
/// Note that the mask name has no relation to names of data items.
void Dataset::setMask(const std::string &maskName, Variable mask) {
  setDims(mask.dims());
  m_masks.insert_or_assign(maskName, std::move(mask));
}

/// Set (insert or replace) an mask for item with given name.
void Dataset::setMask(const std::string &name, const std::string &maskName,
                      Variable mask) {
  scipp::expect::contains(*this, name);
  setDims(mask.dims());
  m_data[name].masks.insert_or_assign(maskName, std::move(mask));
}

void Dataset::setData_impl(const std::string &name, detail::DatasetData &&data,
                           const AttrPolicy attrPolicy) {
  setDims(data.data ? data.data.dims() : data.unaligned->dims);
  const auto replace = contains(name);
  if (replace && attrPolicy == AttrPolicy::Keep)
    data.coords = std::move(m_data[name].coords);
  m_data[name] = std::move(data);
  if (replace)
    rebuildDims();
}

/// Set (insert or replace) data (values, optional variances) with given name.
///
/// Throws if the provided values bring the dataset into an inconsistent state
/// (mismatching dtype, unit, or dimensions). The default is to drop existing
/// attributes, unless AttrPolicy::Keep is specified.
void Dataset::setData(const std::string &name, Variable data,
                      const AttrPolicy attrPolicy) {
  setData_impl(name, detail::DatasetData{std::move(data), {}, {}, {}},
               attrPolicy);
}

/// Private helper for constructor of DataArray and setData
void Dataset::setData(const std::string &name, UnalignedData &&data) {
  setData_impl(
      name,
      detail::DatasetData{
          {}, std::make_unique<UnalignedData>(std::move(data)), {}, {}},
      AttrPolicy::Drop);
}

/// Set (insert or replace) data from a DataArray with a given name, avoiding
/// copies where possible by using std::move.
void Dataset::setData(const std::string &name, DataArray data) {
  // Get the Dataset holder
  auto dataset = DataArray::to_dataset(std::move(data));
  // There can be only one DatasetData item, so get the first one with begin()
  auto item = dataset.m_data.begin();

  for (auto &&[dim, coord] : dataset.m_coords) {
    if (const auto it = m_coords.find(dim); it != m_coords.end())
      core::expect::equals(coord, it->second);
    else
      setCoord(dim, std::move(coord));
  }

  for (auto &&[nm, mask] : dataset.m_masks)
    if (const auto it = m_masks.find(std::string(nm)); it != m_masks.end())
      core::expect::equals(mask, it->second);
    else
      setMask(std::string(nm), std::move(mask));

  if (item->second.data)
    setData(name, std::move(item->second.data));
  else
    setData(name, std::move(*item->second.unaligned));
  for (auto &&[dim, coord] : item->second.coords)
    setCoord(name, dim, std::move(coord));
}

/// Set (insert or replace) data item with given name.
///
/// Coordinates, masks, and attributes of the data array are added to the
/// dataset. Throws if there are existing but mismatching coords, masks, or
/// attributes. Throws if the provided data brings the dataset into an
/// inconsistent state (mismatching dtype, unit, or dimensions).
void Dataset::setData(const std::string &name, const DataArrayConstView &data) {
  if (contains(name) && &m_data[name] == &data.underlying() &&
      data.slices().empty())
    return; // Self-assignment, return early.
  setData(name, DataArray(data));
}

void DataArrayView::setData(Variable data) const {
  if (!slices().empty())
    throw except::SliceError("Cannot set data via slice.");
  if (!m_mutableData->second.data && hasData())
    m_mutableData->second.unaligned->data.setData(std::move(data));
  else
    m_mutableDataset->setData(name(), std::move(data), AttrPolicy::Keep);
}

/// Removes the coordinate for the given dimension.
void Dataset::eraseCoord(const Dim dim) { erase_from_map(m_coords, dim); }

/// Removes unaligned coord with given name from the given item.
void Dataset::eraseCoord(const std::string &name, const Dim dim) {
  scipp::expect::contains(*this, name);
  erase_from_map(m_data[name].coords, dim);
}

/// Remove mask with given name.
void Dataset::eraseMask(const std::string &maskName) {
  erase_from_map(m_masks, maskName);
}

/// Remove mask with given mask name from the given item.
void Dataset::eraseMask(const std::string &name, const std::string &maskName) {
  scipp::expect::contains(*this, name);
  erase_from_map(m_data[name].masks, maskName);
}

/// Return const slice of the dataset along given dimension with given extents.
///
/// This does not make a copy of the data. Instead of view object is returned.
DatasetConstView Dataset::slice(const Slice s) const & {
  return DatasetConstView(*this).slice(s);
}

/// Return slice of the dataset along given dimension with given extents.
///
/// This does not make a copy of the data. Instead of view object is returned.
DatasetView Dataset::slice(const Slice s) & {
  return DatasetView(*this).slice(s);
}

/// Return const slice of the dataset along given dimension with given extents.
///
/// This overload for rvalue reference *this avoids returning a view
/// referencing data that is about to go out of scope and returns a new dataset
/// instead.
Dataset Dataset::slice(const Slice s) const && {
  return Dataset{DatasetConstView(*this).slice(s)};
}

/// Rename dimension `from` to `to`.
void Dataset::rename(const Dim from, const Dim to) {
  if ((from != to) && (m_dims.count(to) != 0))
    throw except::DimensionError("Duplicate dimension.");

  const auto relabel = [from, to](auto &map) {
    auto node = map.extract(from);
    node.key() = to;
    map.insert(std::move(node));
  };
  if (m_dims.count(from) != 0)
    relabel(m_dims);
  if (m_coords.count(from))
    relabel(m_coords);
  for (auto &item : m_coords)
    item.second.rename(from, to);
  for (auto &item : m_masks)
    item.second.rename(from, to);
  for (auto &item : m_data) {
    auto &value = item.second;
    if (value.data)
      value.data.rename(from, to);
    else {
      value.unaligned->dims.relabel(value.unaligned->dims.index(from), to);
      value.unaligned->data.rename(from, to);
    }
    for (auto &coord : value.coords)
      coord.second.rename(from, to);
  }
}

DataArrayConstView::DataArrayConstView(
    const Dataset &dataset, const detail::dataset_item_map::value_type &data,
    const detail::slice_list &slices, VariableView &&view)
    : m_dataset(&dataset), m_data(&data), m_slices(slices) {
  if (view)
    m_view = std::move(view);
  else if (m_data->second.data)
    m_view =
        VariableView(detail::makeSlice(m_data->second.data, this->slices()));
}

/// Return the name of the view.
///
/// The name of the view is equal to the name of the item in a Dataset, or the
/// name of a DataArray. Note that comparison operations ignore the name.
const std::string &DataArrayConstView::name() const noexcept {
  return m_data->first;
}

/// Set the name of a data array.
void DataArray::setName(const std::string &name) {
  auto &map = m_holder.m_data;
  auto node = map.extract(map.begin());
  node.key() = name;
  map.insert(std::move(node));
}

Dimensions DataArrayConstView::parentDims() const noexcept {
  if (underlying().data)
    return underlying().data.dims();
  else if (hasData()) // view of unaligned content
    return underlying().unaligned->data.dims();
  else
    return underlying().unaligned->dims;
}

/// Return an ordered mapping of dimension labels to extents.
Dimensions DataArrayConstView::dims() const noexcept {
  if (hasData())
    return data().dims();
  else {
    Dimensions sliced(m_data->second.unaligned->dims);
    for (const auto &item : slices()) {
      if (item.first.end() == -1)
        sliced.erase(item.first.dim());
      else
        sliced.resize(item.first.dim(), item.first.end() - item.first.begin());
    }
    return sliced;
  }
}

/// Return the dtype of the data.
DType DataArrayConstView::dtype() const {
  return hasData() ? data().dtype()
                   : dataset::unaligned::is_realigned_events(*this)
                         ? event_dtype(unaligned().dtype())
                         : unaligned().dtype();
}

/// Return the unit of the data values.
units::Unit DataArrayConstView::unit() const {
  return hasData() ? data().unit() : unaligned().unit();
}

DataArrayConstView DataArrayConstView::unaligned() const {
  // This needs to combine coords from m_dataset and from the unaligned data. We
  // therefore construct with normal dataset and items pointers and only pass
  // the unaligned data via a variable view.
  if (hasData())
    return DataArrayConstView{};
  auto view = m_data->second.unaligned->data.data();
  detail::do_make_slice(view, slices());
  return {*m_dataset, *m_data, slices(), std::move(view)};
}

DataArrayView DataArrayView::unaligned() const {
  // See also DataArrayConstView::unaligned.
  if (hasData())
    return DataArrayView{};
  auto view = m_mutableData->second.unaligned->data.data();
  detail::do_make_slice(view, slices());
  return {*m_mutableDataset, *m_mutableData, slices(), std::move(view)};
}

/// Set the unit of the data values.
///
/// Throws if there are no data values.
void DataArrayView::setUnit(const units::Unit unit) const {
  if (hasData())
    return data().setUnit(unit);
  throw except::UnalignedError("Realigned data, cannot set unit.");
}

namespace {
auto unaligned_by_dim_slice = [](const auto &item, const Dim dim) {
  const auto &[key, var] = item;
  if constexpr (std::is_same_v<std::decay_t<decltype(key)>, Dim>) {
    const bool is_dimension_coord =
        !contains_events(var) && var.dims().contains(key);
    return !contains_events(var) && var.dims().contains(dim) &&
           (is_dimension_coord ? key == dim : var.dims().inner() == dim);
  } else {
    return false;
  }
};

template <class Items>
void erase_if_unaligned_by_dim_slice(Items &items, const Dim dim) {
  for (auto it = items.begin(); it != items.end();) {
    if (unaligned_by_dim_slice(*it, dim))
      it = items.erase(it);
    else
      ++it;
  }
}

template <class Items, class Slices>
void erase_if_unaligned_by_dim_slices(Items &items, const Slices &slices) {
  for (const auto &slice : slices)
    if (!slice.first.isRange())
      erase_if_unaligned_by_dim_slice(items, slice.first.dim());
}
} // namespace

/*
template <class MapView>
MapView DataArrayConstView::makeView(const bool aligned) const {
  auto map_parent = [](const DataArrayConstView &self) -> auto & {
    if constexpr (std::is_same_v<MapView, CoordsConstView>)
      return self.m_dataset->m_coords;
    else if constexpr (std::is_same_v<MapView, MasksConstView>)
      return self.m_dataset->m_masks;
    else
      return self.m_data->second.attrs; // item attrs, not dataset attrs
  };
  auto items = makeViewItems<MapView>(parentDims(), map_parent(*this));
  if (aligned)
    erase_if_unaligned_by_dim_slices(items, slices());
  if (!m_data->second.data && hasData()) {
    // This is a view of the unaligned content of a realigned data array.
    decltype(*this) unaligned = m_data->second.unaligned->data;
    auto unalignedItems =
        makeViewItems<MapView>(unaligned.dims(), map_parent(unaligned));
    items.insert(unalignedItems.begin(), unalignedItems.end());
  } else {
    // Attributes are per-item in dataset and may exceed data dims, e.g., for
    // storing bounds when creating a copy of a non-range slice of a histogram.
    if constexpr (std::is_same_v<MapView, AttrsConstView>)
      addOrthogonalAttrs(items, parentDims(), map_parent(*this));
  }
  if constexpr (std::is_same_v<MapView, AttrsConstView>)
    addAttrsFromCoords(items, parentDims(), dims(), m_dataset->m_coords);
  if constexpr (std::is_same_v<MapView, MasksConstView>) {
    auto unalignedItems = makeViewItems<MapView>(dims(), m_data->second.masks);
    items.insert(unalignedItems.begin(), unalignedItems.end());
    addMasks(items, parentDims(), dims(), m_dataset->m_masks);
  }
  return MapView(std::move(items), slices());
}

template <class MapView>
MapView DataArrayView::makeView(const bool aligned) const {
  auto map_parent = [](const DataArrayView &self) -> auto & {
    if constexpr (std::is_same_v<MapView, CoordsView>)
      return self.m_mutableDataset->m_coords;
    else if constexpr (std::is_same_v<MapView, MasksView>)
      return self.m_mutableDataset->m_masks;
    else
      return self.m_mutableData->second.attrs; // item attrs, not dataset attrs
  };
  auto items = makeViewItems<MapView>(parentDims(), map_parent(*this));
  if (aligned)
    erase_if_unaligned_by_dim_slices(items, slices());
  const bool is_view_of_unaligned = !m_mutableData->second.data && hasData();
  DataArray *unaligned_ptr{nullptr};
  if (is_view_of_unaligned) {
    unaligned_ptr = &m_mutableData->second.unaligned->data;
    decltype(*this) unaligned = m_mutableData->second.unaligned->data;
    auto unalignedItems =
        makeViewItems<MapView>(unaligned.dims(), map_parent(unaligned));
    items.insert(unalignedItems.begin(), unalignedItems.end());
  } else {
    if constexpr (std::is_same_v<MapView, AttrsView>)
      addOrthogonalAttrs(items, parentDims(), map_parent(*this));
  }
  if constexpr (std::is_same_v<MapView, AttrsView>)
    addAttrsFromCoords(items, parentDims(), dims(), m_mutableDataset->m_coords);
  if constexpr (std::is_same_v<MapView, MasksView>) {
    auto unalignedItems =
        makeViewItems<MapView>(dims(), m_mutableData->second.masks);
    items.insert(unalignedItems.begin(), unalignedItems.end());
    addMasks(items, parentDims(), dims(), m_mutableDataset->m_masks);
  }
  if constexpr (std::is_same_v<MapView, CoordsView>) {
    // Access disabled with nullptr since views of dataset items or slices of
    // data arrays may not set or erase coords.
    return MapView({unaligned_ptr ? m_mutableDataset : nullptr, unaligned_ptr},
                   std::move(items), slices());
  } else {
    // Note: Unlike for CoordAccess this is *not* unconditionally disabled with
    // nullptr since it sets/erase attributes/masks of the *item*.
    return MapView(
        {slices().empty() ? m_mutableDataset : nullptr, &name(), unaligned_ptr},
        std::move(items), slices());
  }
}
*/

CoordsConstView DataArrayConstView::make_coords(const bool aligned) const {
  // Aligned coords (including unaligned by slicing) from dataset
  auto items = makeViewItems(parentDims(), m_dataset->m_coords);
  if (aligned) {
    erase_if_unaligned_by_dim_slices(items, slices());
  } else {
    // Unaligned coords
    // Note no dims passed, include everything
    const auto tmp = makeViewItems(m_data->second.coords);
    items.insert(tmp.begin(), tmp.end());
  }
  if (!m_data->second.data && hasData()) {
    // This is a view of the unaligned content of a realigned data array.
    decltype(*this) unaligned = m_data->second.unaligned->data;
    // TODO handle aligned flag
    const auto tmp =
        makeViewItems(unaligned.parentDims(), unaligned.m_dataset->m_coords);
    items.insert(tmp.begin(), tmp.end());
  }
  return CoordsConstView(std::move(items), slices());
}

/// Return a const view to all coordinates of the data view.
CoordsConstView DataArrayConstView::coords() const noexcept {
  return make_coords(false);
}

/// Return a const view to all aligned coordinates of the data view.
CoordsConstView DataArrayConstView::aligned_coords() const noexcept {
  return make_coords(true);
}

/// Return a const view to all unaligned coordinates of the data view.
CoordsConstView DataArrayConstView::unaligned_coords() const noexcept {
  auto items = coords().items();
  const auto aligned = aligned_coords().items();
  for (auto it = items.begin(); it != items.end();) {
    if (aligned.count(it->first))
      it = items.erase(it);
    else
      ++it;
  }
  return CoordsConstView{std::move(items), slices()};
}

/// Return a const view to all masks of the data view.
MasksConstView DataArrayConstView::masks() const noexcept {
  // Aligned masks from dataset
  auto items = makeViewItems(parentDims(), m_dataset->m_masks);
  // Unaligned masks
  auto tmp = makeViewItems(parentDims(), m_data->second.masks);
  items.insert(tmp.begin(), tmp.end());
  if (!m_data->second.data && hasData()) {
    // This is a view of the unaligned content of a realigned data array.
    decltype(*this) unaligned = m_data->second.unaligned->data;
    tmp = makeViewItems(unaligned.parentDims(), unaligned.m_dataset->m_masks);
    items.insert(tmp.begin(), tmp.end());
  }
  return MasksConstView(std::move(items), slices());
}

/// Return a const view to all coordinates of the data array.
CoordsConstView DataArray::coords() const { return get().coords(); }
/// Return a const view to all aligned coordinates of the data array.
CoordsConstView DataArray::aligned_coords() const {
  return get().aligned_coords();
}
/// Return a const view to all unaligned coordinates of the data array.
CoordsConstView DataArray::unaligned_coords() const {
  return get().unaligned_coords();
}
/// Return a const view to all masks of the data array.
MasksConstView DataArray::masks() const { return get().masks(); }

DataArrayConstView DataArrayConstView::slice(const Slice s) const {
  const auto &dims_ = dims();
  core::expect::validSlice(dims_, s);
  auto tmp(m_slices);
  tmp.emplace_back(s, dims_[s.dim()]);
  if (!m_data->second.data && hasData()) {
    auto view = m_data->second.unaligned->data.data();
    detail::do_make_slice(view, tmp);
    return {*m_dataset, *m_data, std::move(tmp), std::move(view)};
  } else {
    return {*m_dataset, *m_data, std::move(tmp)};
  }
}

DataArrayView::DataArrayView(Dataset &dataset,
                             detail::dataset_item_map::value_type &data,
                             const detail::slice_list &slices,
                             VariableView &&view)
    : DataArrayConstView(dataset, data, slices,
                         data.second.data ? VariableView(detail::makeSlice(
                                                data.second.data, slices))
                                          : std::move(view)),
      m_mutableDataset(&dataset), m_mutableData(&data) {}

DataArrayView DataArrayView::slice(const Slice s) const {
  core::expect::validSlice(dims(), s);
  auto tmp(slices());
  tmp.emplace_back(s, dims()[s.dim()]);
  if (!m_mutableData->second.data && hasData()) {
    auto view = m_mutableData->second.unaligned->data.data();
    detail::do_make_slice(view, tmp);
    return {*m_mutableDataset, *m_mutableData, std::move(tmp), std::move(view)};
  } else {
    return {*m_mutableDataset, *m_mutableData, std::move(tmp)};
  }
}

CoordsView DataArrayView::make_coords(const bool aligned) const {
  // Aligned coords (including unaligned by slicing) from dataset
  auto items = makeViewItems(parentDims(), m_mutableDataset->m_coords);
  if (aligned) {
    erase_if_unaligned_by_dim_slices(items, slices());
  } else {
    // Unaligned coords
    // Note no dims passed, include everything
    const auto tmp = makeViewItems(m_mutableData->second.coords);
    items.insert(tmp.begin(), tmp.end());
  }
  const bool is_view_of_unaligned = !m_mutableData->second.data && hasData();
  DataArray *unaligned_ptr{nullptr};
  if (is_view_of_unaligned) {
    // This is a view of the unaligned content of a realigned data array.
    unaligned_ptr = &m_mutableData->second.unaligned->data;
    decltype(*this) unaligned = m_mutableData->second.unaligned->data;
    // TODO handle aligned flag
    const auto tmp = makeViewItems(unaligned.parentDims(),
                                   unaligned.m_mutableDataset->m_coords);
    items.insert(tmp.begin(), tmp.end());
  }
  return CoordsView(CoordAccess{slices().empty() ? m_mutableDataset : nullptr,
                                &name(), unaligned_ptr},
                    std::move(items), slices());
}

/// Return a view to all coordinates of the data view.
CoordsView DataArrayView::coords() const noexcept { return make_coords(false); }

/// Return a view to all coordinates of the data array.
CoordsView DataArray::coords() {
  auto items = get().coords().items();
  return CoordsView(CoordAccess{&m_holder}, std::move(items), {});
}

/// Return a view to all aligned coordinates of the data view.
CoordsView DataArrayView::aligned_coords() const noexcept {
  return make_coords(true);
}

/// Return a view to all aligned coordinates of the data array.
CoordsView DataArray::aligned_coords() { return m_holder.coords(); }

/// Return a view to all unaligned coordinates of the data view.
CoordsView DataArrayView::unaligned_coords() const noexcept {
  auto items = coords().items();
  const auto aligned = aligned_coords().items();
  for (auto it = items.begin(); it != items.end();) {
    if (aligned.count(it->first))
      it = items.erase(it);
    else
      ++it;
  }
  // TODO need to init last arg of CoordAccess
  return CoordsView{CoordAccess{slices().empty() ? m_mutableDataset : nullptr,
                                &name(), nullptr},
                    std::move(items), slices()};
}

/// Return a view to all unaligned coordinates of the data array.
CoordsView DataArray::unaligned_coords() { return get().unaligned_coords(); }

/// Return a view to all masks of the data view.
MasksView DataArrayView::masks() const noexcept {
  // Aligned masks from dataset
  auto items = makeViewItems(parentDims(), m_mutableDataset->m_masks);
  // Unaligned masks
  auto tmp = makeViewItems(parentDims(), m_mutableData->second.masks);
  items.insert(tmp.begin(), tmp.end());
  const bool is_view_of_unaligned = !m_mutableData->second.data && hasData();
  DataArray *unaligned_ptr{nullptr};
  if (is_view_of_unaligned) {
    // This is a view of the unaligned content of a realigned data array.
    unaligned_ptr = &m_mutableData->second.unaligned->data;
    decltype(*this) unaligned = m_mutableData->second.unaligned->data;
    tmp = makeViewItems(unaligned.parentDims(),
                        unaligned.m_mutableDataset->m_masks);
    items.insert(tmp.begin(), tmp.end());
  }
  return MasksView(MaskAccess{slices().empty() ? m_mutableDataset : nullptr,
                              &name(), unaligned_ptr},
                   std::move(items), slices());
  // return makeView<MasksView>();
}

/// Return a view to all masks of the data array.
MasksView DataArray::masks() {
  // Unlike coords we do not have aligned masks, not using different MaskAccess.
  return get().masks();
}

DataArrayView DataArrayView::assign(const DataArrayConstView &other) const {
  if (&underlying() == &other.underlying() && slices() == other.slices())
    return *this; // Self-assignment, return early.

  expect::coordsAreSuperset(*this, other);
  // TODO here and below: If other has data, we should either fail, or create
  // data.
  if (hasData())
    data().assign(other.data());
  return *this;
}

DataArrayView DataArrayView::assign(const Variable &other) const {
  if (hasData())
    data().assign(other);
  return *this;
}

DataArrayView DataArrayView::assign(const VariableConstView &other) const {
  if (hasData())
    data().assign(other);
  return *this;
}

DatasetView DatasetView::assign(const DatasetConstView &other) const {
  for (const auto &data : other)
    operator[](data.name()).assign(data);
  return *this;
}

DatasetConstView::DatasetConstView(const Dataset &dataset)
    : m_dataset(&dataset) {
  m_items.reserve(dataset.size());
  for (const auto &item : dataset.m_data)
    m_items.emplace_back(DataArrayView(detail::make_item{this}(item)));
}

DatasetView::DatasetView(Dataset &dataset)
    : DatasetConstView(DatasetConstView::makeViewWithEmptyIndexes(dataset)),
      m_mutableDataset(&dataset) {
  m_items.reserve(size());
  for (auto &item : dataset.m_data)
    m_items.emplace_back(detail::make_item{this}(item));
}

/// Return a const view to all coordinates of the dataset slice.
CoordsConstView DatasetConstView::coords() const noexcept {
  auto items = makeViewItems(m_dataset->dimensions(), m_dataset->m_coords);
  erase_if_unaligned_by_dim_slices(items, slices());
  return CoordsConstView(std::move(items), slices());
}

/// Return a view to all coordinates of the dataset slice.
CoordsView DatasetView::coords() const noexcept {
  auto items =
      makeViewItems(m_mutableDataset->dimensions(), m_mutableDataset->m_coords);
  erase_if_unaligned_by_dim_slices(items, slices());
  return CoordsView(CoordAccess(slices().empty() ? m_mutableDataset : nullptr),
                    std::move(items), slices());
}

/// Return a const view to all masks of the dataset slice.
MasksConstView DatasetConstView::masks() const noexcept {
  auto items = makeViewItems(m_dataset->dimensions(), m_dataset->m_masks);
  return MasksConstView(std::move(items), slices());
}

/// Return a view to all masks of the dataset slice.
MasksView DatasetView::masks() const noexcept {
  return MasksView(MaskAccess(slices().empty() ? m_mutableDataset : nullptr),
                   makeViewItems(dimensions(), m_mutableDataset->m_masks),
                   slices());
}

bool DatasetConstView::contains(const std::string &name) const noexcept {
  return find(name) != end();
}

namespace {
template <class T> const auto &getitem(const T &view, const std::string &name) {
  if (auto it = view.find(name); it != view.end())
    return *it;
  throw except::NotFoundError("Expected " + to_string(view) + " to contain " +
                              name + ".");
}
} // namespace

/// Return a const view to data and coordinates with given name.
const DataArrayConstView &
DatasetConstView::operator[](const std::string &name) const {
  return getitem(*this, name);
}

/// Return a view to data and coordinates with given name.
const DataArrayView &DatasetView::operator[](const std::string &name) const {
  return getitem(*this, name);
}

// This is a member so it gets access to a private constructor of DataArrayView.
template <class T>
std::pair<boost::container::small_vector<DataArrayView, 8>, detail::slice_list>
DatasetConstView::slice_items(const T &view, const Slice slice) {
  auto slices = view.slices();
  boost::container::small_vector<DataArrayView, 8> items;
  scipp::index extent = std::numeric_limits<scipp::index>::max();
  for (const auto &item : view) {
    const auto &dims = item.dims();
    if (dims.contains(slice.dim())) {
      items.emplace_back(DataArrayView(item.slice(slice)));
      // In principle data may be on bin edges. The overall dimension is then
      // determined by the extent of data that is *not* on the edges.
      extent = std::min(extent, dims[slice.dim()]);
    }
  }
  if (extent == std::numeric_limits<scipp::index>::max()) {
    // Fallback: Could not determine extent from data (not data that depends on
    // slicing dimension), use `dimensions()` to also consider coords.
    const auto currentDims = view.dimensions();
    core::expect::validSlice(currentDims, slice);
    extent = currentDims.at(slice.dim());
  }
  slices.emplace_back(slice, extent);
  return std::pair{std::move(items), std::move(slices)};
}

/// Return a slice of the dataset view.
///
/// The returned view will not contain references to data items that do not
/// depend on the sliced dimension.
DatasetConstView DatasetConstView::slice(const Slice s) const {
  DatasetConstView sliced;
  sliced.m_dataset = m_dataset;
  std::tie(sliced.m_items, sliced.m_slices) = slice_items(*this, s);
  return sliced;
}

DatasetView DatasetView::slice(const Slice s) const {
  DatasetView sliced;
  sliced.m_dataset = m_dataset;
  sliced.m_mutableDataset = m_mutableDataset;
  std::tie(sliced.m_items, sliced.m_slices) = slice_items(*this, s);
  return sliced;
}

template <class A, class B> bool dataset_equals(const A &a, const B &b) {
  if (a.size() != b.size())
    return false;
  if (a.coords() != b.coords())
    return false;
  if (a.masks() != b.masks())
    return false;
  for (const auto &data : a) {
    try {
      if (data != b[data.name()])
        return false;
    } catch (except::NotFoundError &) {
      return false;
    }
  }
  return true;
}

/// Return true if the datasets have identical content.
bool Dataset::operator==(const Dataset &other) const {
  return dataset_equals(*this, other);
}

/// Return true if the datasets have identical content.
bool Dataset::operator==(const DatasetConstView &other) const {
  return dataset_equals(*this, other);
}

/// Return true if the datasets have identical content.
bool DatasetConstView::operator==(const Dataset &other) const {
  return dataset_equals(*this, other);
}

/// Return true if the datasets have identical content.
bool DatasetConstView::operator==(const DatasetConstView &other) const {
  return dataset_equals(*this, other);
}

/// Return true if the datasets have mismatching content./
bool Dataset::operator!=(const Dataset &other) const {
  return !dataset_equals(*this, other);
}

/// Return true if the datasets have mismatching content.
bool Dataset::operator!=(const DatasetConstView &other) const {
  return !dataset_equals(*this, other);
}

/// Return true if the datasets have mismatching content.
bool DatasetConstView::operator!=(const Dataset &other) const {
  return !dataset_equals(*this, other);
}

/// Return true if the datasets have mismatching content.
bool DatasetConstView::operator!=(const DatasetConstView &other) const {
  return !dataset_equals(*this, other);
}

std::unordered_map<Dim, scipp::index> DatasetConstView::dimensions() const {
  auto base_dims = m_dataset->dimensions();
  // Note current slices are ordered, but NOT unique
  for (const auto &[slice, extents] : m_slices) {
    (void)extents;
    auto it = base_dims.find(slice.dim());
    if (!slice.isRange()) { // For non-range. Erase dimension
      base_dims.erase(it);
    } else {
      it->second = slice.end() -
                   slice.begin(); // Take extent from slice. This is the effect
                                  // that the successful slice range will have
    }
  }
  return base_dims;
}

std::unordered_map<Dim, scipp::index> Dataset::dimensions() const {
  std::unordered_map<Dim, scipp::index> all;
  for (const auto &dim : this->m_dims) {
    all[dim.first] = extents::decodeExtent(dim.second);
  }
  return all;
}

std::map<typename MasksConstView::key_type,
         typename MasksConstView::mapped_type>
union_or(const MasksConstView &currentMasks, const MasksConstView &otherMasks) {
  std::map<typename MasksConstView::key_type,
           typename MasksConstView::mapped_type>
      out;

  for (const auto &[key, item] : currentMasks) {
    out.emplace(key, item);
  }

  for (const auto &[key, item] : otherMasks) {
    const auto it = currentMasks.find(key);
    if (it != currentMasks.end()) {
      out[key] |= item;
    } else {
      out.emplace(key, item);
    }
  }
  return out;
}

void union_or_in_place(const MasksView &currentMasks,
                       const MasksConstView &otherMasks) {
  for (const auto &[key, item] : otherMasks) {
    const auto it = currentMasks.find(key);
    if (it != currentMasks.end()) {
      it->second |= item;
    } else {
      currentMasks.set(key, item);
    }
  }
}
} // namespace scipp::dataset
