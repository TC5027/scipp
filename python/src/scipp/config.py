# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2019 Scipp contributors (https://github.com/scipp)
# @file
# @author Neil Vaytet


class _Plot:

    def __init__(self):

        # The list of default line colors
        self.color = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
                      "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"]

        # The colorbar properties
        self.params = {"cmap": "viridis", "log": False, "vmin": None,
                       "vmax": None, "color": None, "show": True, "cbar": True}

        # The default image height (in pixels)
        self.height = 533

        # The default image width (in pixels)
        self.width = 800

        # Resolution
        self.dpi = 96

        # Aspect ratio for images: "equal" will conserve the aspect ratio,
        # while "auto" will stretch the image to the size of the figure
        self.aspect = "auto"

        # Make list of markers for matplotlib
        self.marker = ["o", "^", "s", "d", "*", "1", "P", "h", "X", "v", "<",
                       ">", "2", "3", "4", "8", "p", "H", "+", "x", "D"]

        # Default line width for 1D plots
        self.linewidth = [1.5]

        # Default line style for 1D non-histogram plots
        self.linestyle = ["none"]

        #
        self.max_plot_size = 100000


class _Colors:

    def __init__(self):

        # The colors for each dataset member used in table and show functions
        self.scheme = {
            "coord": "#dde9af",
            "data": "#ffe680",
            "labels": "#afdde9",
            "attr": "#ff8080",
            "mask": "#cccccc"
        }


plot = _Plot()
colors = _Colors()
