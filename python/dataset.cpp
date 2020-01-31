// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2019 Scipp contributors (https://github.com/scipp)
/// @file
/// @author Simon Heybrock

#include "scipp/core/dataset.h"
#include "scipp/core/except.h"
#include "scipp/core/sort.h"

#include "bind_data_access.h"
#include "bind_operators.h"
#include "bind_slice_methods.h"
#include "detail.h"
#include "pybind11.h"
#include "rename.h"

using namespace scipp;
using namespace scipp::core;

namespace py = pybind11;

/// Helper to provide equivalent of the `items()` method of a Python dict.
template <class T> class items_view {
public:
  items_view(T &obj) : m_obj(&obj) {}
  auto size() const noexcept { return m_obj->size(); }
  auto begin() const { return m_obj->items_begin(); }
  auto end() const { return m_obj->items_end(); }

private:
  T *m_obj;
};
template <class T> items_view(T &)->items_view<T>;

/// Helper to provide equivalent of the `values()` method of a Python dict.
template <class T> class values_view {
public:
  values_view(T &obj) : m_obj(&obj) {}
  auto size() const noexcept { return m_obj->size(); }
  auto begin() const {
    if constexpr (std::is_same_v<typename T::mapped_type, DataArray>)
      return m_obj->begin();
    else
      return m_obj->values_begin();
  }
  auto end() const {
    if constexpr (std::is_same_v<typename T::mapped_type, DataArray>)
      return m_obj->end();
    else
      return m_obj->values_end();
  }

private:
  T *m_obj;
};
template <class T> values_view(T &)->values_view<T>;

/// Helper to provide equivalent of the `keys()` method of a Python dict.
template <class T> class keys_view {
public:
  keys_view(T &obj) : m_obj(&obj) {}
  auto size() const noexcept { return m_obj->size(); }
  auto begin() const { return m_obj->keys_begin(); }
  auto end() const { return m_obj->keys_end(); }

private:
  T *m_obj;
};
template <class T> keys_view(T &)->keys_view<T>;

template <template <class> class View, class T>
void bind_helper_view(py::module &m, const std::string &name) {
  std::string suffix;
  if (std::is_same_v<View<T>, items_view<T>>)
    suffix = "_items_view";
  if (std::is_same_v<View<T>, values_view<T>>)
    suffix = "_values_view";
  if (std::is_same_v<View<T>, keys_view<T>>)
    suffix = "_keys_view";
  py::class_<View<T>>(m, (name + suffix).c_str())
      .def(py::init([](T &obj) { return View{obj}; }))
      .def("__len__", &View<T>::size)
      .def("__iter__",
           [](View<T> &self) {
             return py::make_iterator(self.begin(), self.end(),
                                      py::return_value_policy::move);
           },
           py::return_value_policy::move, py::keep_alive<0, 1>());
}

template <class T, class ConstT>
void bind_mutable_proxy(py::module &m, const std::string &name) {
  py::class_<ConstT>(m, (name + "ConstProxy").c_str());
  py::class_<T, ConstT> proxy(m, (name + "Proxy").c_str());
  proxy.def("__len__", &T::size)
      .def("__getitem__", &T::operator[], py::return_value_policy::move,
           py::keep_alive<0, 1>())
      .def("__setitem__",
           [](T &self, const typename T::key_type key,
              const VariableConstProxy &var) { self.set(key, var); })
      // This additional setitem allows us to do things like
      // d.attrs["a"] = scipp.detail.move(scipp.Variable())
      .def("__setitem__",
           [](T &self, const typename T::key_type key, MoveableVariable &mvar) {
             self.set(key, std::move(mvar.var));
           })
      .def("__delitem__", &T::erase, py::call_guard<py::gil_scoped_release>())
      .def("__iter__",
           [](T &self) {
             return py::make_iterator(self.keys_begin(), self.keys_end(),
                                      py::return_value_policy::move);
           },
           py::keep_alive<0, 1>())
      .def("keys", [](T &self) { return keys_view(self); },
           py::keep_alive<0, 1>(), R"(view on self's keys)")
      .def("values", [](T &self) { return values_view(self); },
           py::keep_alive<0, 1>(), R"(view on self's values)")
      .def("items", [](T &self) { return items_view(self); },
           py::return_value_policy::move, py::keep_alive<0, 1>(),
           R"(view on self's items)")
      .def("__contains__", &T::contains);
  bind_comparison<T>(proxy);
}

template <class T, class... Ignored>
void bind_coord_properties(py::class_<T, Ignored...> &c) {
  // For some reason the return value policy and/or keep-alive policy do not
  // work unless we wrap things in py::cpp_function.
  c.def_property_readonly(
      "coords",
      py::cpp_function([](T &self) { return self.coords(); },
                       py::return_value_policy::move, py::keep_alive<0, 1>()),
      R"(
      Dict of coordinates.)");
  c.def_property_readonly(
      "labels",
      py::cpp_function([](T &self) { return self.labels(); },
                       py::return_value_policy::move, py::keep_alive<0, 1>()),
      R"(
      Dict of labels.

      Labels are very similar to coordinates, except that they are identified
      using custom names instead of dimension labels.)");
  c.def_property_readonly("masks",
                          py::cpp_function([](T &self) { return self.masks(); },
                                           py::return_value_policy::move,
                                           py::keep_alive<0, 1>()),
                          R"(
      Dict of masks.)");
  c.def_property_readonly("attrs",
                          py::cpp_function([](T &self) { return self.attrs(); },
                                           py::return_value_policy::move,
                                           py::keep_alive<0, 1>()),
                          R"(
      Dict of attributes.)");
}

template <class T, class... Ignored>
void bind_dataset_proxy_methods(py::class_<T, Ignored...> &c) {
  c.def("__len__", &T::size);
  c.def("__repr__", [](const T &self) { return to_string(self); });
  c.def("__iter__",
        [](const T &self) {
          return py::make_iterator(self.keys_begin(), self.keys_end(),
                                   py::return_value_policy::move);
        },
        py::return_value_policy::move, py::keep_alive<0, 1>());
  c.def("keys", [](T &self) { return keys_view(self); },
        py::return_value_policy::move, py::keep_alive<0, 1>(),
        R"(view on self's keys)");
  c.def("values", [](T &self) { return values_view(self); },
        py::return_value_policy::move, py::keep_alive<0, 1>(),
        R"(view on self's values)");
  c.def("items", [](T &self) { return items_view(self); },
        py::return_value_policy::move, py::keep_alive<0, 1>(),
        R"(view on self's items)");
  c.def("__getitem__",
        [](T &self, const std::string &name) { return self[name]; },
        py::keep_alive<0, 1>());
  c.def("__contains__", &T::contains);
  c.def("copy", [](const T &self) { return Dataset(self); },
        py::call_guard<py::gil_scoped_release>(), "Return a (deep) copy.");
  c.def("__copy__", [](const T &self) { return Dataset(self); },
        py::call_guard<py::gil_scoped_release>(), "Return a (deep) copy.");
  c.def("__deepcopy__",
        [](const T &self, const py::dict &) { return Dataset(self); },
        py::call_guard<py::gil_scoped_release>(), "Return a (deep) copy.");
  c.def_property_readonly("dims",
                          [](const T &self) {
                            std::vector<Dim> dims;
                            for (const auto &dim : self.dimensions()) {
                              dims.push_back(dim.first);
                            }
                            return dims;
                          },
                          R"(List of dimensions.)",
                          py::return_value_policy::move);
  c.def_property_readonly("shape",
                          [](const T &self) {
                            auto shape = py::list();
                            for (const auto &dim : self.dimensions()) {
                              if (dim.second == Dimensions::Sparse) {
                                shape.append(py::none());
                              } else {
                                shape.append(dim.second);
                              }
                            }
                            return shape;
                          },
                          R"(List of shapes.)", py::return_value_policy::move);
}

template <class T, class... Ignored>
void bind_data_array_properties(py::class_<T, Ignored...> &c) {
  c.def_property_readonly("name", &T::name, R"(The name of the held data.)");
  c.def("__repr__", [](const T &self) { return to_string(self); });
  c.def("copy", [](const T &self) { return DataArray(self); },
        py::call_guard<py::gil_scoped_release>(), "Return a (deep) copy.");
  c.def("__copy__", [](const T &self) { return DataArray(self); },
        py::call_guard<py::gil_scoped_release>(), "Return a (deep) copy.");
  c.def("__deepcopy__",
        [](const T &self, const py::dict &) { return DataArray(self); },
        py::call_guard<py::gil_scoped_release>(), "Return a (deep) copy.");
  c.def_property(
      "data",
      py::cpp_function(
          [](T &self) {
            return self.hasData() ? py::cast(self.data()) : py::none();
          },
          py::return_value_policy::move, py::keep_alive<0, 1>()),
      [](T &self, const VariableConstProxy &data) { self.data().assign(data); },
      R"(Underlying data item.)");
  bind_coord_properties(c);
  bind_comparison<DataConstProxy>(c);
  bind_data_properties(c);
  bind_slice_methods(c);
  bind_in_place_binary<DataProxy>(c);
  bind_in_place_binary<VariableConstProxy>(c);
  bind_binary<Dataset>(c);
  bind_binary<DatasetProxy>(c);
  bind_binary<DataProxy>(c);
  bind_binary<VariableConstProxy>(c);
}

void init_dataset(py::module &m) {
  py::class_<Slice>(m, "Slice");

  bind_helper_view<items_view, Dataset>(m, "Dataset");
  bind_helper_view<items_view, DatasetProxy>(m, "DatasetProxy");
  bind_helper_view<items_view, CoordsProxy>(m, "CoordsProxy");
  bind_helper_view<items_view, LabelsProxy>(m, "LabelsProxy");
  bind_helper_view<items_view, MasksProxy>(m, "MasksProxy");
  bind_helper_view<items_view, AttrsProxy>(m, "AttrsProxy");
  bind_helper_view<keys_view, Dataset>(m, "Dataset");
  bind_helper_view<keys_view, DatasetProxy>(m, "DatasetProxy");
  bind_helper_view<keys_view, CoordsProxy>(m, "CoordsProxy");
  bind_helper_view<keys_view, LabelsProxy>(m, "LabelsProxy");
  bind_helper_view<keys_view, MasksProxy>(m, "MasksProxy");
  bind_helper_view<keys_view, AttrsProxy>(m, "AttrsProxy");
  bind_helper_view<values_view, Dataset>(m, "Dataset");
  bind_helper_view<values_view, DatasetProxy>(m, "DatasetProxy");
  bind_helper_view<values_view, CoordsProxy>(m, "CoordsProxy");
  bind_helper_view<values_view, LabelsProxy>(m, "LabelsProxy");
  bind_helper_view<values_view, MasksProxy>(m, "MasksProxy");
  bind_helper_view<values_view, AttrsProxy>(m, "AttrsProxy");

  bind_mutable_proxy<CoordsProxy, CoordsConstProxy>(m, "Coords");
  bind_mutable_proxy<LabelsProxy, LabelsConstProxy>(m, "Labels");
  bind_mutable_proxy<MasksProxy, MasksConstProxy>(m, "Masks");
  bind_mutable_proxy<AttrsProxy, AttrsConstProxy>(m, "Attrs");

  py::class_<DataArray> dataArray(m, "DataArray", R"(
    Named variable with associated coords, labels, and attributes.)");
  dataArray.def(py::init<const DataConstProxy &>());
  dataArray.def(
      py::init<std::optional<Variable>, std::map<Dim, Variable>,
               std::map<std::string, Variable>, std::map<std::string, Variable>,
               std::map<std::string, Variable>>(),
      py::arg("data") = std::nullopt,
      py::arg("coords") = std::map<Dim, Variable>{},
      py::arg("labels") = std::map<std::string, Variable>{},
      py::arg("masks") = std::map<std::string, Variable>{},
      py::arg("attrs") = std::map<std::string, Variable>{});

  py::class_<DataConstProxy>(m, "DataConstProxy")
      .def(py::init<const DataArray &>());

  py::class_<DataProxy, DataConstProxy> dataProxy(m, "DataProxy", R"(
        Proxy for DataArray, representing a sliced view onto a DataArray, or an item of a Dataset;
        Mostly equivalent to DataArray, see there for details.)");
  dataProxy.def(py::init<DataArray &>());

  bind_data_array_properties(dataArray);
  bind_data_array_properties(dataProxy);

  py::class_<DatasetConstProxy>(m, "DatasetConstProxy")
      .def(py::init<const Dataset &>());
  py::class_<DatasetProxy, DatasetConstProxy> datasetProxy(m, "DatasetProxy",
                                                           R"(
        Proxy for Dataset, representing a sliced view onto a Dataset;
        Mostly equivalent to Dataset, see there for details.)");
  datasetProxy.def(py::init<Dataset &>());

  py::class_<Dataset> dataset(m, "Dataset", R"(
    Dict of data arrays with aligned dimensions.)");

  dataset.def(py::init<const std::map<std::string, DataConstProxy> &>())
      .def(py::init<const DataConstProxy &>())
      .def(py::init([](const std::map<std::string, VariableConstProxy> &data,
                       const std::map<Dim, VariableConstProxy> &coords,
                       const std::map<std::string, VariableConstProxy> &labels,
                       const std::map<std::string, VariableConstProxy> &masks,
                       const std::map<std::string, VariableConstProxy> &attrs) {
             return Dataset(data, coords, labels, masks, attrs);
           }),
           py::arg("data") = std::map<std::string, VariableConstProxy>{},
           py::arg("coords") = std::map<Dim, VariableConstProxy>{},
           py::arg("labels") = std::map<std::string, VariableConstProxy>{},
           py::arg("masks") = std::map<std::string, VariableConstProxy>{},
           py::arg("attrs") = std::map<std::string, VariableConstProxy>{})
      .def(py::init([](const DatasetProxy &other) { return Dataset{other}; }))
      .def("__setitem__",
           [](Dataset &self, const std::string &name,
              const VariableConstProxy &data) { self.setData(name, data); })
      .def("__setitem__",
           [](Dataset &self, const std::string &name, MoveableVariable &mvar) {
             self.setData(name, std::move(mvar.var));
           })
      .def("__setitem__",
           [](Dataset &self, const std::string &name,
              const DataConstProxy &data) { self.setData(name, data); })
      .def("__setitem__",
           [](Dataset &self, const std::string &name, MoveableDataArray &mdat) {
             self.setData(name, std::move(mdat.data));
           })
      .def("__delitem__", &Dataset::erase,
           py::call_guard<py::gil_scoped_release>())
      .def(
          "clear", &Dataset::clear,
          R"(Removes all data (preserving coordinates, attributes, labels and masks.).)");
  datasetProxy.def("__setitem__",
                   [](const DatasetProxy &self, const std::string &name,
                      const DataConstProxy &data) { self[name].assign(data); });

  bind_dataset_proxy_methods(dataset);
  bind_dataset_proxy_methods(datasetProxy);

  bind_coord_properties(dataset);
  bind_coord_properties(datasetProxy);

  bind_slice_methods(dataset);
  bind_slice_methods(datasetProxy);

  bind_comparison<Dataset>(dataset);
  bind_comparison<DatasetProxy>(dataset);
  bind_comparison<Dataset>(datasetProxy);
  bind_comparison<DatasetProxy>(datasetProxy);

  bind_in_place_binary<Dataset>(dataset);
  bind_in_place_binary<DatasetProxy>(dataset);
  bind_in_place_binary<DataProxy>(dataset);
  bind_in_place_binary<VariableConstProxy>(dataset);
  bind_in_place_binary<Dataset>(datasetProxy);
  bind_in_place_binary<DatasetProxy>(datasetProxy);
  bind_in_place_binary<VariableConstProxy>(datasetProxy);
  bind_in_place_binary<DataProxy>(datasetProxy);
  bind_in_place_binary_scalars(dataset);
  bind_in_place_binary_scalars(datasetProxy);
  bind_in_place_binary_scalars(dataArray);
  bind_in_place_binary_scalars(dataProxy);

  bind_binary<Dataset>(dataset);
  bind_binary<DatasetProxy>(dataset);
  bind_binary<DataProxy>(dataset);
  bind_binary<VariableConstProxy>(dataset);
  bind_binary<Dataset>(datasetProxy);
  bind_binary<DatasetProxy>(datasetProxy);
  bind_binary<DataProxy>(datasetProxy);
  bind_binary<VariableConstProxy>(datasetProxy);

  dataArray.def("rename_dims", &rename_dims<DataArray>, py::arg("dims_dict"),
                "Rename dimensions.");
  dataset.def("rename_dims", &rename_dims<Dataset>, py::arg("dims_dict"),
              "Rename dimensions.");

  m.def("concatenate",
        py::overload_cast<const DataConstProxy &, const DataConstProxy &,
                          const Dim>(&concatenate),
        py::arg("x"), py::arg("y"), py::arg("dim"),
        py::call_guard<py::gil_scoped_release>(), R"(
        Concatenate input data array along the given dimension.

        Concatenates the data, coords, labels and masks of the data array.
        Coords, labels and masks for any but the given dimension are required to match and are copied to the output without changes.

        :param x: First DataArray.
        :param y: Second DataArray.
        :param dim: Dimension along which to concatenate.
        :raises: If the dtype or unit does not match, or if the dimensions and shapes are incompatible.
        :return: New data array containing all data, coords, labels, and masks of the input arrays.
        :rtype: DataArray)");

  m.def("concatenate",
        py::overload_cast<const DatasetConstProxy &, const DatasetConstProxy &,
                          const Dim>(&concatenate),
        py::arg("x"), py::arg("y"), py::arg("dim"),
        py::call_guard<py::gil_scoped_release>(), R"(
        Concatenate input datasets along the given dimension.

        Concatenate all cooresponding items in the input datasets.
        The output contains only items that are present in both inputs.

        :param x: First Dataset.
        :param y: Second Dataset.
        :param dim: Dimension along which to concatenate.
        :raises: If the dtype or unit does not match, or if the dimensions and shapes are incompatible.
        :return: New dataset.
        :rtype: Dataset)");

  m.def("histogram",
        [](const DataConstProxy &ds, const Variable &bins) {
          return core::histogram(ds, bins);
        },
        py::arg("x"), py::arg("bins"), py::call_guard<py::gil_scoped_release>(),
        R"(Returns a new DataArray with values in bins for sparse dims.

        :param x: Data to histogram.
        :param bins: Bin edges.
        :return: Histogramed data.
        :rtype: DataArray)");

  m.def("histogram",
        [](const DataConstProxy &ds, const VariableConstProxy &bins) {
          return core::histogram(ds, bins);
        },
        py::arg("x"), py::arg("bins"), py::call_guard<py::gil_scoped_release>(),
        R"(Returns a new DataArray with values in bins for sparse dims.

        :param x: Data to histogram.
        :param bins: Bin edges.
        :return: Histogramed data.
        :rtype: DataArray)");

  m.def("histogram",
        [](const Dataset &ds, const VariableConstProxy &bins) {
          return core::histogram(ds, bins);
        },
        py::arg("x"), py::arg("bins"), py::call_guard<py::gil_scoped_release>(),
        R"(Returns a new Dataset with values in bins for sparse dims.

        :param x: Data to histogram.
        :param bins: Bin edges.
        :return: Histogramed data.
        :rtype: Dataset)");

  m.def("histogram",
        [](const Dataset &ds, const Variable &bins) {
          return core::histogram(ds, bins);
        },
        py::arg("x"), py::arg("bins"), py::call_guard<py::gil_scoped_release>(),
        R"(Returns a new Dataset with values in bins for sparse dims.

        :param x: Data to histogram.
        :param bins: Bin edges.
        :return: Histogramed data.
        :rtype: Dataset)");

  m.def("merge",
        [](const DatasetConstProxy &lhs, const DatasetConstProxy &rhs) {
          return core::merge(lhs, rhs);
        },
        py::arg("lhs"), py::arg("rhs"),
        py::call_guard<py::gil_scoped_release>(), R"(
        Union of two datasets.

        :param lhs: First Dataset.
        :param rhs: Second Dataset.
        :raises: If there are conflicting items with different content.
        :return: A new dataset that contains the union of all data items, coords, labels, masks and attributes.
        :rtype: Dataset)");

  m.def("sum", py::overload_cast<const DataConstProxy &, const Dim>(&sum),
        py::arg("x"), py::arg("dim"), py::call_guard<py::gil_scoped_release>(),
        R"(
        Element-wise sum over the specified dimension.

        :param x: Data to sum.
        :param dim: Dimension over which to sum.
        :raises: If the dimension does not exist, or if the dtype cannot be summed, e.g., if it is a string
        :seealso: :py:class:`scipp.mean`
        :return: New data array containing the sum.
        :rtype: DataArray)");

  m.def("sum", py::overload_cast<const DatasetConstProxy &, const Dim>(&sum),
        py::arg("x"), py::arg("dim"), py::call_guard<py::gil_scoped_release>(),
        R"(
        Element-wise sum over the specified dimension.

        :param x: Data to sum.
        :param dim: Dimension over which to sum.
        :raises: If the dimension does not exist, or if the dtype cannot be summed, e.g., if it is a string
        :seealso: :py:class:`scipp.mean`
        :return: New dataset containing the sum for each data item.
        :rtype: Dataset)");

  m.def("mean", py::overload_cast<const DataConstProxy &, const Dim>(&mean),
        py::arg("x"), py::arg("dim"), py::call_guard<py::gil_scoped_release>(),
        R"(
        Element-wise mean over the specified dimension, if variances are present, the new variance is computated as standard-deviation of the mean.

        See the documentation for the mean of a Variable for details in the computation of the ouput variance.

        :param x: Data to calculate mean of.
        :param dim: Dimension over which to calculate mean.
        :raises: If the dimension does not exist, or if the dtype cannot be summed, e.g., if it is a string
        :seealso: :py:class:`scipp.mean`
        :return: New data array containing the mean for each data item.
        :rtype: DataArray)");

  m.def("mean", py::overload_cast<const DatasetConstProxy &, const Dim>(&mean),
        py::arg("x"), py::arg("dim"), py::call_guard<py::gil_scoped_release>(),
        R"(
        Element-wise mean over the specified dimension, if variances are present, the new variance is computated as standard-deviation of the mean.

        See the documentation for the mean of a Variable for details in the computation of the ouput variance.

        :param x: Data to calculate mean of.
        :param dim: Dimension over which to calculate mean.
        :raises: If the dimension does not exist, or if the dtype cannot be summed, e.g., if it is a string
        :seealso: :py:class:`scipp.mean`
        :return: New dataset containing the mean for each data item.
        :rtype: Dataset)");

  m.def("rebin",
        py::overload_cast<const DataConstProxy &, const Dim,
                          const VariableConstProxy &>(&rebin),
        py::arg("x"), py::arg("dim"), py::arg("bins"),
        py::call_guard<py::gil_scoped_release>(), R"(
        Rebin a dimension of a data array.

        :param x: Data to rebin.
        :param dim: Dimension to rebin over.
        :param bins: New bin edges.
        :raises: If data cannot be rebinned, e.g., if the unit is not counts, or the existing coordinate is not a bin-edge coordinate.
        :return: A new data array with data rebinned according to the new coordinate.
        :rtype: DataArray)");

  m.def("rebin",
        py::overload_cast<const DatasetConstProxy &, const Dim,
                          const VariableConstProxy &>(&rebin),
        py::arg("x"), py::arg("dim"), py::arg("bins"),
        py::call_guard<py::gil_scoped_release>(), R"(
        Rebin a dimension of a dataset.

        :param x: Data to rebin.
        :param dim: Dimension to rebin over.
        :param bins: New bin edges.
        :raises: If data cannot be rebinned, e.g., if the unit is not counts, or the existing coordinate is not a bin-edge coordinate.
        :return: A new dataset with data rebinned according to the new coordinate.
        :rtype: Dataset)");

  m.def("sort",
        py::overload_cast<const DataConstProxy &, const VariableConstProxy &>(
            &sort),
        py::arg("x"), py::arg("key"), py::call_guard<py::gil_scoped_release>(),
        R"(Sort data array along a dimension by a sort key.

        :raises: If the key is invalid, e.g., if it has not exactly one dimension, or if its dtype is not sortable.
        :return: New sorted data array.
        :rtype: DataArray)");

  m.def(
      "sort", py::overload_cast<const DataConstProxy &, const Dim &>(&sort),
      py::arg("x"), py::arg("key"), py::call_guard<py::gil_scoped_release>(),
      R"(Sort data array along a dimension by the coordinate values for that dimension.

      :raises: If the key is invalid, e.g., if it has not exactly one dimension, or if its dtype is not sortable.
      :return: New sorted data array.
      :rtype: DataArray)");

  m.def(
      "sort",
      py::overload_cast<const DataConstProxy &, const std::string &>(&sort),
      py::arg("x"), py::arg("key"), py::call_guard<py::gil_scoped_release>(),
      R"(Sort data array along a dimension by the label values for the given key.

      :raises: If the key is invalid, e.g., if it has not exactly one dimension, or if its dtype is not sortable.
      :return: New sorted data array.
      :rtype: DataArray)");

  m.def(
      "sort",
      py::overload_cast<const DatasetConstProxy &, const VariableConstProxy &>(
          &sort),
      py::arg("x"), py::arg("key"), py::call_guard<py::gil_scoped_release>(),
      R"(Sort dataset along a dimension by a sort key.

        :raises: If the key is invalid, e.g., if it has not exactly one dimension, or if its dtype is not sortable.
        :return: New sorted dataset.
        :rtype: Dataset)");

  m.def(
      "sort", py::overload_cast<const DatasetConstProxy &, const Dim &>(&sort),
      py::arg("x"), py::arg("key"), py::call_guard<py::gil_scoped_release>(),
      R"(Sort dataset along a dimension by the coordinate values for that dimension.

      :raises: If the key is invalid, e.g., if it has not exactly one dimension, or if its dtype is not sortable.
      :return: New sorted dataset.
      :rtype: Dataset)");

  m.def(
      "sort",
      py::overload_cast<const DatasetConstProxy &, const std::string &>(&sort),
      py::arg("x"), py::arg("key"), py::call_guard<py::gil_scoped_release>(),
      R"(Sort dataset along a dimension by the label values for the given key.

      :raises: If the key is invalid, e.g., if it has not exactly one dimension, or if its dtype is not sortable.
      :return: New sorted dataset.
      :rtype: Dataset)");

  m.def("combine_masks",
        [](const MasksConstProxy &msk, const std::vector<Dim> &labels,
           const std::vector<scipp::index> &shape) {
          return core::masks_merge_if_contained(
              msk, core::Dimensions(labels, shape));
        },
        py::call_guard<py::gil_scoped_release>(), R"(
        Combine all masks into a single one following the OR operation.
        This requires a masks proxy as an input, followed by the dimension
        labels and shape of the Variable/DataArray. The labels and the shape
        are used to create a Dimensions object. The function then iterates
        through the masks proxy and combines only the masks that have all
        their dimensions contained in the Variable/DataArray Dimensions.

        :return: A new variable that contains the union of all masks.
        :rtype: Variable)");

  m.def("reciprocal",
        [](const DataConstProxy &self) { return reciprocal(self); },
        py::arg("x"), py::call_guard<py::gil_scoped_release>(), R"(
        Element-wise reciprocal.

        :return: Reciprocal of the input values.
        :rtype: DataArray)");

  py::implicitly_convertible<DataArray, DataConstProxy>();
  py::implicitly_convertible<DataArray, DataProxy>();
  py::implicitly_convertible<Dataset, DatasetConstProxy>();
}
