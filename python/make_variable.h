// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2020 Scipp contributors (https://github.com/scipp)
/// @file
/// @author Simon Heybrock

#include "scipp/units/unit.h"

// #include "scipp/core/dtype.h"
// #include "scipp/core/except.h"
// #include "scipp/core/tag_util.h"

// #include "scipp/variable/comparison.h"
// #include "scipp/variable/operations.h"
// #include "scipp/variable/transform.h"
#include "scipp/variable/variable.h"

// #include "scipp/dataset/dataset.h"
// #include "scipp/dataset/sort.h"

// #include "bind_data_access.h"
// #include "bind_operators.h"
// #include "bind_slice_methods.h"
// #include "dtype.h"
// #include "numpy.h"
// #include "py_object.h"
#include "pybind11.h"
// #include "rename.h"

using namespace scipp;
using namespace scipp::variable;

namespace py = pybind11;

template <class T> struct MakeVariable {
  static Variable apply(const std::vector<Dim> &labels, py::array values,
                        const std::optional<py::array> &variances,
                        const units::Unit unit) {
    // Pybind11 converts py::array to py::array_t for us, with all sorts of
    // automatic conversions such as integer to double, if required.
    py::array_t<T> valuesT(values);
    py::buffer_info info = valuesT.request();
    Dimensions dims(labels, {info.shape.begin(), info.shape.end()});
    auto var = variances
                   ? makeVariable<T>(Dimensions{dims}, Values{}, Variances{})
                   : makeVariable<T>(Dimensions(dims));
    copy_flattened<T>(valuesT, var.template values<T>());
    if (variances) {
      py::array_t<T> variancesT(*variances);
      info = variancesT.request();
      core::expect::equals(
          dims, Dimensions(labels, {info.shape.begin(), info.shape.end()}));
      copy_flattened<T>(variancesT, var.template variances<T>());
    }
    var.setUnit(unit);
    return var;
  }
};

template <class T> struct MakeVariableDefaultInit {
  static Variable apply(const std::vector<Dim> &labels,
                        const std::vector<scipp::index> &shape,
                        const units::Unit unit, const bool variances) {
    Dimensions dims(labels, shape);
    auto var = variances
                   ? makeVariable<T>(Dimensions{dims}, Values{}, Variances{})
                   : makeVariable<T>(Dimensions{dims});
    var.setUnit(unit);
    return var;
  }
};

template <class ST> struct MakeODFromNativePythonTypes {
  template <class T> struct Maker {
    static Variable apply(const units::Unit unit, const ST &value,
                          const std::optional<ST> &variance) {
      auto var = variance ? makeVariable<T>(Values{T(value)},
                                            Variances{T(variance.value())})
                          : makeVariable<T>(Values{T(value)});
      var.setUnit(unit);
      return var;
    }
  };

  static Variable make(const units::Unit unit, const ST &value,
                       const std::optional<ST> &variance,
                       const py::object &dtype) {
    return core::CallDType<double, float, int64_t, int32_t, bool>::apply<Maker>(
        scipp_dtype(dtype), unit, value, variance);
  }
};

template <class T>
Variable init_1D_no_variance(const std::vector<Dim> &labels,
                             const std::vector<scipp::index> &shape,
                             const std::vector<T> &values,
                             const units::Unit &unit) {
  Variable var;
  var = makeVariable<T>(Dims(labels), Shape(shape),
                        Values(values.begin(), values.end()));
  var.setUnit(unit);
  return var;
}

template <class T>
auto do_init_0D(const T &value, const std::optional<T> &variance,
                const units::Unit &unit) {
  using Elem = std::conditional_t<std::is_same_v<T, py::object>,
                                  scipp::python::PyObject, T>;
  Variable var;
  if (variance)
    var = makeVariable<Elem>(Values{value}, Variances{*variance});
  else
    var = makeVariable<Elem>(Values{value});
  var.setUnit(unit);
  return var;
}

Variable doMakeVariable(const std::vector<Dim> &labels, py::array &values,
                        std::optional<py::array> &variances,
                        const units::Unit unit, const py::object &dtype) {
  // Use custom dtype, otherwise dtype of data.
  const auto dtypeTag =
      dtype.is_none() ? scipp_dtype(values.dtype()) : scipp_dtype(dtype);

  if (labels.size() == 1 && !variances) {
    if (dtypeTag == core::dtype<std::string>) {
      std::vector<scipp::index> shape(values.shape(),
                                      values.shape() + values.ndim());
      return init_1D_no_variance(labels, shape,
                                 values.cast<std::vector<std::string>>(), unit);
    }

    if (dtypeTag == core::dtype<Eigen::Vector3d> ||
        dtypeTag == core::dtype<Eigen::Quaterniond>) {
      std::vector<scipp::index> shape(values.shape(),
                                      values.shape() + values.ndim() - 1);
      if (dtypeTag == core::dtype<Eigen::Vector3d>) {
        return init_1D_no_variance(
            labels, shape, values.cast<std::vector<Eigen::Vector3d>>(), unit);
      } else {
        const auto &arr = values.cast<std::vector<std::vector<double>>>();
        std::vector<Eigen::Quaterniond> qvec;
        qvec.reserve(arr.size());
        for (size_t i = 0; i < arr.size(); i++)
          qvec.emplace_back(arr[i].data());
        return init_1D_no_variance(labels, shape, qvec, unit);
      }
    }
  }

  return core::CallDType<double, float, int64_t, int32_t,
                         bool>::apply<MakeVariable>(dtypeTag, labels, values,
                                                    variances, unit);
}

Variable makeVariableDefaultInit(const std::vector<Dim> &labels,
                                 const std::vector<scipp::index> &shape,
                                 const units::Unit unit, py::object &dtype,
                                 const bool variances) {
  return core::CallDType<
      double, float, int64_t, int32_t, bool, event_list<double>,
      event_list<float>, event_list<int64_t>, event_list<int32_t>, DataArray,
      Dataset, Eigen::Vector3d,
      Eigen::Quaterniond>::apply<MakeVariableDefaultInit>(scipp_dtype(dtype),
                                                          labels, shape, unit,
                                                          variances);
}
