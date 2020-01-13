# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2019 Scipp contributors (https://github.com/scipp)
# @author Neil Vaytet

from ..config import plot as config
from .tools import parse_params
from ..utils import name_with_unit, value_to_string
from .._scipp.core import combine_masks, Variable, Dim, dtype


# Other imports
import numpy as np
import matplotlib.ticker as ticker
import math


class Slicer:

    def __init__(self, scipp_obj_dict=None, data_array=None, axes=None,
                 values=None, variances=None, masks=None, cmap=None, log=None,
                 vmin=None, vmax=None, color=None, button_options=None,
                 volume=False, aspect=None, max_plot_size=512):

        import ipywidgets as widgets

        self.scipp_obj_dict = scipp_obj_dict
        if self.scipp_obj_dict is None:
            self.scipp_obj_dict = {data_array.name: data_array}
        self.data_array = data_array

        if max_plot_size is None:
            self.max_plot_size = config.max_plot_size
        else:
            self.max_plot_size = max_plot_size

        # Member container for dict output
        self.members = dict(widgets=dict(sliders=dict(), togglebuttons=dict(),
                            togglebutton=dict(), buttons=dict(),
                            labels=dict()))

        # Parse parameters for values, variances and masks
        self.params = dict()
        globs = {"cmap": cmap, "log": log, "vmin": vmin, "vmax": vmax,
                 "color": color}

        self.params["values"] = parse_params(params=values, globs=globs,
                                             array=self.data_array.values)

        self.params["variances"] = {"show": False}
        if self.data_array.variances is not None:
            self.params["variances"].update(
                parse_params(params=variances, defaults={"show": False},
                             globs=globs,
                             array=np.sqrt(self.data_array.variances)))

        self.params["masks"] = parse_params(
            params=masks, defaults={"cmap": "gray", "cbar": False},
            globs=globs)
        self.params["masks"]["show"] = (self.params["masks"]["show"] and
                                        len(self.data_array.masks) > 0)
        if self.params["masks"]["show"]:
            self.masks = combine_masks(self.data_array.masks,
                                       self.data_array.dims,
                                       self.data_array.shape)

        self.labels = self.data_array.labels
        self.shapes = dict(zip(self.data_array.dims, self.data_array.shape))
        # Save aspect ratio setting
        self.aspect = aspect
        if self.aspect is None:
            self.aspect = config.aspect

        # Size of the slider coordinate arrays
        self.slider_nx = dict()
        # Store coordinates of dimensions that will be in sliders
        self.slider_x = dict()
        # Store ticklabels for a dimension
        self.slider_ticks = dict()
        # Store labels for sliders if any
        self.slider_labels = dict()

        # Process axes dimensions
        if axes is None:
            axes = self.data_array.dims
        # Protect against duplicate entries in axes
        if len(axes) != len(set(axes)):
            raise RuntimeError("Duplicate entry in axes: {}".format(axes))
        # Iterate through axes and collect dimensions
        for ax in axes:
            dim, lab, var, ticks = self.axis_label_and_ticks(ax)
            if (lab is not None) and (dim in axes):
                raise RuntimeError("The dimension of the labels cannot also "
                                   "be specified as another axis.")
            self.slider_labels[dim] = lab
            self.slider_x[dim] = var
            self.slider_ticks[dim] = ticks
            self.slider_nx[dim] = self.shapes[dim]
        self.ndim = len(self.slider_nx)

        # Save information on histograms
        self.histograms = dict()
        for name, var in self.scipp_obj_dict.items():
            self.histograms[name] = dict()
            for dim, x in self.slider_x.items():
                indx = var.dims.index(dim)
                self.histograms[name][dim] = var.shape[indx] == x.shape[0] - 1

        # Initialise list for VBox container
        self.vbox = []

        # Initialise slider and label containers
        self.lab = dict()
        self.slider = dict()
        self.buttons = dict()
        self.showhide = dict()
        self.button_axis_to_dim = dict()
        # Default starting index for slider
        indx = 0

        # Now begin loop to construct sliders
        button_values = [None] * (self.ndim - len(button_options)) + \
            button_options[::-1]
        for i, dim in enumerate(self.slider_x.keys()):
            key = str(dim)
            # If this is a 3d projection, place slices half-way
            if len(button_options) == 3 and (not volume):
                indx = (self.slider_nx[dim] - 1) // 2
            if self.slider_labels[dim] is not None:
                descr = self.slider_labels[dim]
            else:
                descr = str(dim)
            # Add an IntSlider to slide along the z dimension of the array
            self.slider[dim] = widgets.IntSlider(
                value=indx,
                min=0,
                max=self.slider_nx[dim] - 1,
                step=1,
                description=descr,
                continuous_update=True,
                readout=False,
                disabled=((i >= self.ndim-len(button_options)) and
                          ((len(button_options) < 3) or volume)))
            labvalue = self.make_slider_label(self.slider_x[dim], indx)
            if self.ndim == len(button_options):
                self.slider[dim].layout.display = 'none'
                labvalue = descr
            # Add a label widget to display the value of the z coordinate
            self.lab[dim] = widgets.Label(value=labvalue)
            # Add one set of buttons per dimension
            self.buttons[dim] = widgets.ToggleButtons(
                options=button_options, description='',
                value=button_values[i],
                disabled=False,
                button_style='',
                style={"button_width": "70px"})
            if button_values[i] is None:
                button_style = ""
            else:
                button_style = "success"
                self.button_axis_to_dim[button_values[i].lower()] = dim
            setattr(self.buttons[dim], "dim", dim)
            setattr(self.buttons[dim], "old_value", self.buttons[dim].value)
            setattr(self.slider[dim], "dim", dim)

            if self.ndim == 1:
                self.buttons[dim].layout.display = 'none'
                self.lab[dim].layout.display = 'none'

            if (len(button_options) == 3) and (not volume):
                self.showhide[dim] = widgets.Button(
                    description="hide",
                    disabled=(button_values[i] is None),
                    button_style=button_style,
                    layout={'width': "70px"}
                )
                setattr(self.showhide[dim], "dim", dim)
                setattr(self.showhide[dim], "value",
                        button_values[i] is not None)
                # Add observer to show/hide buttons
                self.showhide[dim].on_click(self.update_showhide)
                self.members["widgets"]["buttons"][key] = self.showhide[dim]

            # Add observer to buttons
            self.buttons[dim].on_msg(self.update_buttons)
            # Add an observer to the slider
            self.slider[dim].observe(self.update_slice, names="value")
            # Add the row of slider + buttons
            row = [self.slider[dim], self.lab[dim], self.buttons[dim]]
            if (len(button_options) == 3) and (not volume):
                row += [widgets.HTML(value="&nbsp;&nbsp;&nbsp;&nbsp;"),
                        self.showhide[dim]]
            self.vbox.append(widgets.HBox(row))

            # Construct members object
            self.members["widgets"]["sliders"][key] = self.slider[dim]
            self.members["widgets"]["togglebuttons"][key] = self.buttons[dim]
            self.members["widgets"]["labels"][key] = self.lab[dim]

        if self.params["masks"]["show"]:
            self.masks_button = widgets.ToggleButton(
                value=self.params["masks"]["show"],
                description="Hide masks" if self.params["masks"]["show"] else
                            "Show masks",
                disabled=False, button_style="")
            self.masks_button.observe(self.toggle_masks, names="value")
            self.vbox += [self.masks_button]
            self.members["widgets"]["togglebutton"]["masks"] = \
                self.masks_button



        return

    def make_slider_label(self, var, indx):
        return name_with_unit(var=var, name=value_to_string(var.values[indx]))

    def mask_to_float(self, mask, var):
        return np.where(mask, var, None).astype(np.float)

    def axis_label_and_ticks(self, axis):
        """
        Get dimensions and label (if present) from requested axis
        """
        ticks = None
        if isinstance(axis, Dim):
            dim = axis
            lab = None
            make_fake_coord = False
            fake_unit = None
            if not self.data_array.coords.__contains__(dim):
                make_fake_coord = True
            else:
                tp = self.data_array.coords[dim].dtype
                if tp == dtype.vector_3_float64:
                    make_fake_coord = True
                    ticks = {"formatter": lambda x: "(" + ",".join(
                                 [value_to_string(item, precision=2)
                                  for item in x]) + ")"}
                elif tp == dtype.string:
                    make_fake_coord = True
                    ticks = {"formatter": lambda x: x}
                if make_fake_coord:
                    fake_unit = self.data_array.coords[dim].unit
                    ticks["coord"] = self.data_array.coords[dim]
            if make_fake_coord:
                args = {"values": np.arange(self.shapes[dim])}
                if fake_unit is not None:
                    args["unit"] = fake_unit
                var = Variable([dim], **args)
            else:
                var = self.data_array.coords[dim]
        elif isinstance(axis, str):
            # By convention, the last dim of the labels is the inner dimension,
            # but note that for now two-dimensional labels are not supported in
            # the plotting.
            dim = self.data_array.labels[axis].dims[-1]
            lab = axis
            var = self.data_array.labels[lab]
        else:
            raise RuntimeError("Unsupported axis found in 'axes': {}. This "
                               "must be either a Scipp dimension "
                               "or a string.".format(axis))
        return dim, lab, var, ticks

    def get_custom_ticks(self, ax, dim, xy="x"):
        """
        Return a list of string to be used as axis tick labels in the case of
        strings or vectors as one axis coordinate.
        """
        getticks = getattr(ax, "get_{}ticks".format(xy))
        xticks = getticks()
        if xticks[2] - xticks[1] < 1:
            ax.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
            xticks = getticks()
        new_ticks = [""] * len(xticks)
        for i, x in enumerate(xticks):
            if x >= 0 and x < self.slider_nx[dim]:
                new_ticks[i] = self.slider_ticks[dim]["formatter"](
                    self.slider_ticks[dim]["coord"].values[int(x)])
        return new_ticks

    def make_zoom_layers(self):
        # Construct layers for zooming
        data_size = 1
        ndims = 0
        self.nlayers = np.Inf
        for dim, sl in self.slider.items():
            if sl.disabled:
                self.nlayers = min(self.nlayers, int(math.log(self.shapes[dim] / self.max_plot_size) / math.log(2.0)))
        #     if self.shapes[dim]
        #     if sl.disabled:
        #         data_size *= self.shapes[dim]
        #         ndims += 1
        # nlayers = min(10, (data_size / config.max_plot_size)**(1.0/ndims))
        print("nlayers", self.nlayers)

