// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2020 Scipp contributors (https://github.com/scipp)
/// @file
/// @author Simon Heybrock
#include "docstring.h"
#include "pybind11.h"

#include "scipp/variable/misc_operations.h"
#include "scipp/variable/variable.h"

using namespace scipp;
using namespace scipp::variable;

namespace py = pybind11;

template <class T> void bind_reshape(py::module &m) {
  auto doc = Docstring()
          .description("Reshape a variable.")
          .raises("If the volume of the old shape is not equal to the volume "
                  "of the new shape.")
          .returns("New variable with requested dimension labels and shape.")
          .rtype<T>()
          .template param<T>("x", "Input to reshape.")
          .param("dims", "List of new dimensions.", "list")
          .param("shape", "New extents in each dimension.", "tuple");
  m.def(
      "reshape",
      [](const typename T::const_view_type &x, const std::vector<Dim> &labels,
         const py::tuple &shape) {
        Dimensions dims(labels, shape.cast<std::vector<scipp::index>>());
        return x.reshape(dims);
      },
      py::arg("x"), py::arg("dims"), py::arg("shape"),
      py::call_guard<py::gil_scoped_release>(),
      doc.c_str());
  m.def(
      "reshape",
      [](const typename T::const_view_type &x, const std::vector<Dim> &labels,
         const py::list &shape) {
        Dimensions dims(labels, shape.cast<std::vector<scipp::index>>());
        return x.reshape(dims);
      },
      py::arg("x"), py::arg("dims"), py::arg("shape"),
      py::call_guard<py::gil_scoped_release>(),
      doc.param("shape", "New extents in each dimension.", "list").c_str());
}

template <class T> void bind_split(py::module &m) {
  m.def("split",
      [](const typename T::const_view_type &x, const Dim dim,
         const std::vector<scipp::index> &indices){
        return split(x, dim, indices);
       },
        py::arg("x"), py::arg("dim"), py::arg("indices"),
        py::call_guard<py::gil_scoped_release>(),
        Docstring()
          .description("Split a Variable along a given Dimension.")
          .returns("The input data, split along the given dimension.")
          .rtype<T>()
          .template param<T>("x", "Input to split.")
          .param("dim", "Dimension to split.", "str")
          .param("indices", "Indices where to split the input.", "list").c_str());
}


template <class T> void bind_transpose(py::module &m) {
  // auto doc = Docstring()
  //         .description("Transpose a variable.")
  //         // .raises("If the volume of the old shape is not equal to the volume "
  //         //         "of the new shape.")
  //         .returns("A copy of the input with requested dimension labels.")
  //         .rtype<T>()
  //         .template param<T>("x", "Input to transpose.")
  //         .param("dims", "New dimension order.", "list");
  m.def(
      "transpose",
      [](const typename T::const_view_type &x, const std::vector<Dim> &labels) {
        // Note transpose returns a ConstView, so we make a copy.
        return copy(x.transpose(labels));
      },
      py::arg("x"), py::arg("dims"),
      py::call_guard<py::gil_scoped_release>(),
      py::keep_alive<0, 1>(),
      Docstring()
          .description("Transpose a data structure.")
          .raises("If the requested dimension labels are not found in the input.")
          .returns("A copy of the input with new order for dimension labels.")
          .rtype<T>()
          .template param<T>("x", "Input to transpose.")
          .param("dims", "New dimension order.", "list").c_str());
}

void init_shape(py::module &m) {
  bind_reshape<Variable>(m);

  bind_split<Variable>(m);

  bind_transpose<Variable>(m);
}
