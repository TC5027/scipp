# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2020 Scipp contributors (https://github.com/scipp)
# @author Neil Vaytet

# Scipp imports
from .. import config
from .render import render_plot
from .slicer import Slicer
from .tools import centers_to_edges, edges_to_centers, parse_params
from ..utils import name_with_unit
from .._scipp import core as sc

# Other imports
import numpy as np
import ipywidgets as widgets
import matplotlib.pyplot as plt
import warnings
from scipy.interpolate import griddata
from scipy.stats import binned_statistic_2d


def plot_2d(scipp_obj_dict=None,
            axes=None,
            values=None,
            variances=None,
            masks=None,
            filename=None,
            figsize=None,
            mpl_axes=None,
            aspect=None,
            cmap=None,
            log=False,
            vmin=None,
            vmax=None,
            color=None,
            logx=False,
            logy=False,
            logxy=False):
    """
    Plot a 2D slice through a N dimensional dataset. For every dimension above
    2, a slider is created to adjust the position of the slice in that
    particular dimension.
    """

    sv = Slicer2d(scipp_obj_dict=scipp_obj_dict,
                  axes=axes,
                  values=values,
                  variances=variances,
                  masks=masks,
                  mpl_axes=mpl_axes,
                  aspect=aspect,
                  cmap=cmap,
                  log=log,
                  vmin=vmin,
                  vmax=vmax,
                  color=color,
                  logx=logx or logxy,
                  logy=logy or logxy)

    if mpl_axes is None:
        render_plot(figure=sv.fig, widgets=sv.vbox, filename=filename)

    return sv.members


class Slicer2d(Slicer):
    def __init__(self,
                 scipp_obj_dict=None,
                 axes=None,
                 values=None,
                 variances=None,
                 masks=None,
                 mpl_axes=None,
                 aspect=None,
                 cmap=None,
                 log=None,
                 vmin=None,
                 vmax=None,
                 color=None,
                 logx=False,
                 logy=False):

        super().__init__(scipp_obj_dict=scipp_obj_dict,
                         axes=axes,
                         values=values,
                         variances=variances,
                         masks=masks,
                         cmap=cmap,
                         log=log,
                         vmin=vmin,
                         vmax=vmax,
                         color=color,
                         aspect=aspect,
                         button_options=['X', 'Y'])

        self.members.update({"images": {}, "colorbars": {}})
        self.extent = {"x": [1, 2], "y": [1, 2]}
        self.logx = logx
        self.logy = logy
        self.vminmax = {"vmin": vmin, "vmax": vmax}
        self.global_vmin = np.Inf
        self.global_vmax = np.NINF
        self.image_resolution = 0.64 * config.plot.dpi / 96.0
        self.image_resolution = [int(self.image_resolution * config.plot.width),
        int(self.image_resolution * config.plot.height)]

        # Get or create matplotlib axes
        self.fig = None
        cax = [None] * (1 + self.params["variances"][self.name]["show"])
        if mpl_axes is not None:
            if isinstance(mpl_axes, dict):
                ax = [None, None]
                for key, val in mpl_axes.items():
                    if key == "ax" or key == "ax_values":
                        ax[0] = val
                    if key == "cax" or key == "cax_values":
                        cax[0] = val
                    if key == "ax_variances":
                        ax[1] = val
                    if key == "cax_variances":
                        cax[1] = val
            else:
                # Case where only a single axis is given
                ax = [mpl_axes]
        else:
            self.fig, ax = plt.subplots(
                1,
                1 + self.params["variances"][self.name]["show"],
                figsize=(config.plot.width / config.plot.dpi,
                         config.plot.height /
                         (1.0 + self.params["variances"][self.name]["show"]) /
                         config.plot.dpi),
                dpi=config.plot.dpi,
                sharex=True,
                sharey=True)
            if not self.params["variances"][self.name]["show"]:
                ax = [ax]

        self.ax = dict()
        self.cax = dict()
        self.im = dict()
        self.cbar = dict()

        self.ax["values"] = ax[0]
        self.cax["values"] = cax[0]
        panels = ["values"]
        if self.params["variances"][self.name]["show"]:
            self.ax["variances"] = ax[1]
            self.cax["variances"] = cax[1]
            panels.append("variances")

        extent_array = np.array(list(self.extent.values())).flatten()
        for key in panels:
            if self.params[key][self.name]["show"]:
                self.im[key] = self.ax[key].imshow(
                    [[1.0, 1.0], [1.0, 1.0]],
                    norm=self.params[key][self.name]["norm"],
                    extent=extent_array,
                    origin="lower",
                    aspect=self.aspect,
                    interpolation="nearest",
                    cmap=self.params[key][self.name]["cmap"])
                self.ax[key].set_title(self.name if key ==
                                       "values" else "std dev.")
                if self.params[key][self.name]["cbar"]:
                    self.cbar[key] = plt.colorbar(self.im[key],
                                                  ax=self.ax[key],
                                                  cax=self.cax[key])
                    self.cbar[key].ax.set_ylabel(
                        name_with_unit(var=self.data_array, name=""))
                if self.cax[key] is None:
                    self.cbar[key].ax.yaxis.set_label_coords(-1.1, 0.5)
                self.members["images"][key] = self.im[key]
                self.members["colorbars"][key] = self.cbar[key]
                if self.params["masks"][self.name]["show"]:
                    self.im[self.get_mask_key(key)] = self.ax[key].imshow(
                        [[1.0, 1.0], [1.0, 1.0]],
                        extent=extent_array,
                        norm=self.params[key][self.name]["norm"],
                        origin="lower",
                        interpolation="nearest",
                        aspect=self.aspect,
                        cmap=self.params["masks"][self.name]["cmap"])
                if self.logx:
                    self.ax[key].set_xscale("log")
                if self.logy:
                    self.ax[key].set_yscale("log")

        # Call update_slice once to make the initial image
        self.update_axes()
        self.update_slice(None)
        self.vbox = widgets.VBox(self.vbox)
        self.vbox.layout.align_items = 'center'
        self.members["fig"] = self.fig
        self.members["ax"] = self.ax

        return

    def update_buttons(self, owner, event, dummy):
        toggle_slider = False
        if not self.slider[owner.dim].disabled:
            toggle_slider = True
            self.slider[owner.dim].disabled = True
        for dim, button in self.buttons.items():
            if (button.value == owner.value) and (dim != owner.dim):
                if self.slider[dim].disabled:
                    button.value = owner.old_value
                else:
                    button.value = None
                button.old_value = button.value
                if toggle_slider:
                    self.slider[dim].disabled = False
        owner.old_value = owner.value
        self.update_axes()
        self.update_slice(None)

        return

    def update_axes(self):
        # Go through the buttons and select the right coordinates for the axes
        axparams = {"x": {}, "y": {}}
        for dim, button in self.buttons.items():
            if self.slider[dim].disabled:
                but_val = button.value.lower()
                # xc = self.slider_x[self.name][dim].values
                if not self.histograms[self.name][dim]:
                    xc = self.slider_x[self.name][dim].values
                    if self.slider_nx[self.name][dim] < 2:
                        dx = 0.5 * abs(xc[0])
                        if dx == 0.0:
                            dx = 0.5
                        xmin = xc[0] - dx
                        xmax = xc[0] + dx
                        axparams[but_val]["xmin"] = xmin
                        axparams[but_val]["xmax"] = xmax
                    else:
                        xmin = 1.5 * xc[0] - 0.5 * xc[1]
                        xmax = 1.5 * xc[-1] - 0.5 * xc[-2]
                    self.extent[but_val] = [xmin, xmax]
                else:
                    self.extent[but_val] = self.slider_x[
                        self.name][dim].values[[0, -1]].astype(np.float)

                axparams[but_val]["lims"] = self.extent[but_val].copy()
                if getattr(self,
                           "log" + but_val) and (self.extent[but_val][0] <= 0):
                    if not self.histograms[self.name][dim]:
                        new_x = centers_to_edges(xc)
                    else:
                        new_x = edges_to_centers(
                            self.slider_x[self.name][dim].values)
                    axparams[but_val]["lims"][0] = new_x[np.searchsorted(
                        new_x, 0)]
                axparams[but_val]["labels"] = name_with_unit(
                    self.slider_x[self.name][dim], name=str(dim))
                axparams[but_val]["dim"] = dim

        extent_array = np.array(list(self.extent.values())).flatten()

        # Image re-sampling
        # res = 128
        # nx = int(self.image_resolution * config.plot.width)
        # ny = int(self.image_resolution * config.plot.height)
        # print(nx, ny)
        # # xmin = 0.0
        # # xmax = x[-1]
        # # dx = (extent_array[1] - extent_array[0])/float(nx)
        # self.img_xe = np.linspace(extent_array[0], extent_array[1], nx+1)
        # self.img_xc = edges_to_centers(self.img_xe)
        # self.img_ye = np.linspace(extent_array[2], extent_array[3], ny+1)
        # self.img_yc = edges_to_centers(self.img_ye)
        # # ymin = 0.0
        # # ymax = y[-1]
        # # dy = (ymax - ymin)/float(ny)
        # # ye = np.linspace(ymin, ymax, ny+1)
        # # yc = np.linspace(ymin + 0.5*dy, ymax - 0.5*dy, ny)
        # # self.img_xg, self.img_yg = np.meshgrid(self.img_xc, self.img_yc)
        # self.img_xg, self.img_yg = np.meshgrid(self.img_xc, self.img_yc)

        # xx, yy = np.meshgrid(self.slider_x[self.name][axparams['x']["dim"]].values,
        #                      self.slider_x[self.name][axparams['y']["dim"]].values)
        #                      # vslice.coords[button_dims[0]].values)
        # # print(button_dims)
        # self.xflat = xx.ravel()
        # self.yflat = yy.ravel()
        # self.xbinwidth = sc.Variable(dims=[axparams['x']["dim"]], values=np.ediff1d(self.slider_x[self.name][axparams['x']["dim"]].values))
        # self.ybinwidth = sc.Variable(dims=[axparams['y']["dim"]], values=np.ediff1d(self.slider_x[self.name][axparams['y']["dim"]].values))

        self.xrebin = sc.Variable(dims=[axparams['x']["dim"]], values=np.linspace(extent_array[0], extent_array[1], self.image_resolution[0]+1), unit=self.slider_x[self.name][axparams['x']["dim"]].unit)
        self.yrebin = sc.Variable(dims=[axparams['y']["dim"]], values=np.linspace(extent_array[2], extent_array[3], self.image_resolution[1]+1), unit=self.slider_x[self.name][axparams['y']["dim"]].unit)
        
        # self.slider_x[self.name][axparams['x']["dim"]].values
        

        for key in self.ax.keys():
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", category=UserWarning)
                self.im[key].set_extent(extent_array)
                if self.params["masks"][self.name]["show"]:
                    self.im[self.get_mask_key(key)].set_extent(extent_array)
                self.ax[key].set_xlim(axparams["x"]["lims"])
                self.ax[key].set_ylim(axparams["y"]["lims"])
            self.ax[key].set_xlabel(axparams["x"]["labels"])
            self.ax[key].set_ylabel(axparams["y"]["labels"])
            for xy, param in axparams.items():
                if self.slider_ticks[self.name][param["dim"]] is not None:
                    getattr(self.ax[key], "set_{}ticklabels".format(xy))(
                        self.get_custom_ticks(ax=self.ax[key],
                                              dim=param["dim"],
                                              xy=xy))
        return

    def update_slice(self, change):
        """
        Slice data according to new slider value.
        """
        vslice = self.data_array
        if self.params["masks"][self.name]["show"]:
            mslice = self.masks
        # Slice along dimensions with active sliders
        button_dims = [None, None]
        for dim, val in self.slider.items():
            if not val.disabled:
                self.lab[dim].value = self.make_slider_label(
                    self.slider_x[self.name][dim], val.value)
                vslice = vslice[val.dim, val.value]
                # At this point, after masks were combined, all their
                # dimensions should be contained in the data_array.dims.
                if self.params["masks"][self.name]["show"]:
                    mslice = mslice[val.dim, val.value]
            else:
                button_dims[self.buttons[dim].value.lower() == "x"] = val.dim

        # Check if dimensions of arrays agree, if not, plot the transpose
        slice_dims = vslice.dims
        transp = slice_dims != button_dims

        if self.params["masks"][self.name]["show"]:
            shape_list = [self.shapes[self.name][bdim] for bdim in button_dims]
            # Use scipp's automatic broadcast functionality to broadcast
            # lower dimension masks to higher dimensions.
            # TODO: creating a Variable here could become expensive when
            # sliders are being used. We could consider performing the
            # automatic broadcasting once and store it in the Slicer class,
            # but this could create a large memory overhead if the data is
            # large.
            # Here, the data is at most 2D, so having the Variable creation
            # and broadcasting should remain cheap.
            msk = sc.Variable(dims=button_dims,
                           values=np.ones(shape_list, dtype=np.int32))
            msk *= sc.Variable(dims=mslice.dims,
                            values=mslice.values.astype(np.int32))

        autoscale_cbar = False
        if vslice.unaligned is not None:
            vslice = sc.histogram(vslice)
            autoscale_cbar = True

        # if not self.histograms[self.name][dim]:
        #             xc = self.slider_x[self.name][dim].values

        # xx, yy = np.meshgrid(vslice.coords[button_dims[1]].values,
        #                      vslice.coords[button_dims[0]].values)
        # # print(button_dims)
        # xflat = xx.ravel()
        # yflat = yy.ravel()

        # to_process = {"values": vslice.values.ravel()}
        # if "variances" in self.ax.keys():
        #     to_process["variances"] = np.sqrt(vslice.variances.ravel())
        # # print(self.img_ye)
        # # print(self.img_xe)

        # results, y_edges, x_edges, bin_number = binned_statistic_2d(
        #     x=self.yflat, y=self.xflat, values=list(to_process.values()), statistic='mean', bins=[self.img_ye, self.img_xe])

        # subset = np.where(np.isfinite(results[0].ravel()))
        # points = np.transpose([self.img_xg.ravel()[subset], self.img_yg.ravel()[subset]])

        # vslice =  sc.rebin(vslice * self.xbinwidth * self.ybinwidth, self.xrebin.dims[0], self.xrebin)
        # vslice =  sc.rebin(vslice * self.ybinwidth, self.yrebin.dims[0], self.yrebin)
        # vslice =  sc.histogram(vslice, self.xrebin)

        # print(vslice.variances)
        if not self.histograms[self.name][self.xrebin.dims[0]] or not self.histograms[self.name][self.yrebin.dims[0]]:
            vslice = vslice.copy()
        if not self.histograms[self.name][self.xrebin.dims[0]]:
            xe = centers_to_edges(vslice.coords[self.xrebin.dims[0]].values)
            # # sc.reshape(vslice.coords[self.xrebin.dims[0]], self.xrebin.dims, [vslice.coords[self.xrebin.dims[0]].shape[0] + 1 ])
            # print(vslice.coords[self.xrebin.dims[0]].shape)
            # # sc.reshape(vslice.coords[self.xrebin.dims[0]], ['x'], [vslice.coords[self.xrebin.dims[0]].shape[0] + 1 ])
            # print(vslice.coords[self.xrebin.dims[0]].shape[0] + 1)
            # print(len(xe), vslice.coords[self.xrebin.dims[0]].shape, len(vslice.coords[self.xrebin.dims[0]].values))
            # vslice = vslice.copy()
            vslice.coords[self.xrebin.dims[0]] = sc.Variable(dims=self.xrebin.dims, values=xe, unit=vslice.coords[self.xrebin.dims[0]].unit)
        if not self.histograms[self.name][self.yrebin.dims[0]]:
            ye = centers_to_edges(vslice.coords[self.yrebin.dims[0]].values)
            # vslice = vslice.copy()
            vslice.coords[self.yrebin.dims[0]] = sc.Variable(dims=self.yrebin.dims, values=ye, unit=vslice.coords[self.yrebin.dims[0]].unit)
            # sc.reshape(vslice.coords[self.yrebin.dims[0]], self.yrebin.dims, [vslice.coords[self.yrebin.dims[0]].shape[0] + 1])
            # vslice.coords[self.yrebin.dims[0]].values = ye
        vslice =  sc.resample(vslice, self.xrebin.dims[0], self.xrebin, "max")
        vslice =  sc.resample(vslice, self.yrebin.dims[0], self.yrebin, "max")


        for i, key in enumerate(self.ax.keys()):
            arr = getattr(vslice, key)
            # arr = griddata(points, results[i].ravel()[subset], (self.img_xg, self.img_yg), method='nearest')
            # arr = results[i]
            if key == "variances":
                arr = np.sqrt(arr)
            if transp:
                arr = arr.T
            self.im[key].set_data(arr)
            if autoscale_cbar:
                cbar_params = parse_params(globs=self.vminmax,
                                           array=arr,
                                           min_val=self.global_vmin,
                                           max_val=self.global_vmax)
                self.global_vmin = cbar_params["vmin"]
                self.global_vmax = cbar_params["vmax"]
                self.params[key][self.name]["norm"] = cbar_params["norm"]
                self.im[key].set_norm(self.params[key][self.name]["norm"])
            if self.params["masks"][self.name]["show"]:
                self.im[self.get_mask_key(key)].set_data(
                    self.mask_to_float(msk.values, arr))
                self.im[self.get_mask_key(key)].set_norm(
                    self.params[key][self.name]["norm"])

        return

    def toggle_masks(self, change):
        for key in self.ax.keys():
            self.im[key + "_masks"].set_visible(change["new"])
        change["owner"].description = "Hide masks" if change["new"] else \
            "Show masks"
        return

    def get_mask_key(self, key):
        return key + "_masks"
