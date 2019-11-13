# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2019 Scipp contributors (https://github.com/scipp)
# @author Neil Vaytet

# Scipp imports
from ..config import plot as config
from .render import render_plot
from .tools import edges_to_centers, get_1d_axes

# Other imports
import numpy as np
import plotly.graph_objs as go


def plot_1d(input_data, backend=None, logx=False, logy=False, logxy=False,
            color=None, filename=None, axes=None, show_masks=True):
    """
    Plot a 1D spectrum.

    Input is a dictionary containing a list of DataProxy.
    If the coordinate of the x-axis contains bin edges, then a bar plot is
    made.

    TODO: find a more general way of handling arguments to be sent to plotly,
    probably via a dictionay of arguments
    """

    data = []
    ymax = 1.0e-30
    for i, (name, var) in enumerate(sorted(input_data)):

        xlab, ylab, x, y = get_1d_axes(var, axes, name)
        ymax = max(ymax, np.amax(y))

        nx = x.shape[0]
        ny = y.shape[0]
        histogram = False
        if nx == ny + 1:
            histogram = True

        # Define trace
        trace = dict(x=x, y=y, name=name, type="scattergl")
        if histogram:
            trace["line"] = {"shape": "hv"}
            trace["y"] = np.concatenate((trace["y"], [0.0]))
            trace["fill"] = "tozeroy"
            trace["mode"] = "lines"
        if color is not None:
            trace["marker"] = {"color": color[i]}
        # Include variance if present
        if var.variances is not None:
            err_dict = dict(
                    type="data",
                    array=np.sqrt(var.variances),
                    visible=True,
                    color=color[i])
            if histogram:
                trace2 = dict(x=edges_to_centers(x), y=y, showlegend=False,
                              type="scattergl", mode="markers",
                              error_y=err_dict,
                              marker={"color": color[i]})
                data.append(trace2)
            else:
                trace["error_y"] = err_dict

        data.append(trace)

    layout = dict(
        xaxis=dict(title=xlab),
        yaxis=dict(title=ylab),
        showlegend=True,
        legend=dict(x=0.0, y=1.15, orientation="h"),
        height=config.height
    )
    if histogram:
        layout["barmode"] = "overlay"
    if logx or logxy:
        layout["xaxis"]["type"] = "log"
    if logy or logxy:
        layout["yaxis"]["type"] = "log"

    # Add masks
    for name, var in sorted(input_data.masks):
        trace = dict(x=x, y=var.values * ymax,
                              type="scattergl", mode="none",
                              fill="tozeroy")
        data.append(trace)


    fig = go.Figure(data=data, layout=layout)
    render_plot(static_fig=fig, interactive_fig=fig, backend=backend,
                filename=filename)
    return
