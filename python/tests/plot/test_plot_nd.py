# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2020 Scipp contributors (https://github.com/scipp)
# @file
# @author Neil Vaytet

import numpy as np
import pytest
import scipp as sc
from plot_helper import make_dense_dataset, make_events_dataset
from scipp.plot import plot

# Prevent figure from being displayed when running the tests
import matplotlib.pyplot as plt
plt.ioff()

# TODO: For now we are just checking that the plot does not throw any errors.
# In the future it would be nice to check the output by either comparing
# checksums or by using tools like squish.


def test_plot_sliceviewer():
    d = make_dense_dataset(ndim=3)
    plot(d)


def test_plot_sliceviewer_with_two_sliders():
    d = make_dense_dataset(ndim=4)
    plot(d)


def test_plot_sliceviewer_with_axes():
    d = make_dense_dataset(ndim=3)
    plot(d, axes=['x', 'tof', 'y'])


def test_plot_sliceviewer_with_labels():
    d = make_dense_dataset(ndim=3, labels=True)
    plot(d, axes=['x', 'y', "somelabels"])


def test_plot_sliceviewer_with_binedges():
    d = make_dense_dataset(ndim=3, binedges=True)
    plot(d)


def test_plot_projection_3d():
    d = make_dense_dataset(ndim=3)
    plot(d, projection="3d")


def test_plot_projection_3d_with_labels():
    d = make_dense_dataset(ndim=3, labels=True)
    plot(d, projection="3d", axes=['x', 'y', "somelabels"])


def test_plot_projection_3d_with_bin_edges():
    d = make_dense_dataset(ndim=3, binedges=True)
    plot(d, projection="3d")


def test_plot_projection_3d_with_masks():
    d = make_dense_dataset(ndim=3, masks=True)
    plot(d, projection="3d")


def test_plot_projection_3d_with_vectors():
    N = 1000
    M = 100
    theta = np.random.random(N) * np.pi
    phi = np.random.random(N) * 2.0 * np.pi
    r = 10.0 + (np.random.random(N) - 0.5)
    x = r * np.sin(theta) * np.sin(phi)
    y = r * np.sin(theta) * np.cos(phi)
    z = r * np.cos(theta)
    tof = np.arange(M).astype(np.float)
    a = np.arange(M * N).reshape([M, N]) * np.sin(y)
    d = sc.Dataset()
    d.coords['xyz'] = sc.Variable(['xyz'],
                                  values=np.array([x, y, z]).T,
                                  dtype=sc.dtype.vector_3_float64)
    d.coords['tof'] = sc.Variable(['tof'], values=tof)
    d['a'] = sc.Variable(['tof', 'xyz'], values=a)
    plot(d, projection="3d", positions="xyz")


def test_plot_convenience_methods():
    d = make_dense_dataset(ndim=3)
    sc.plot.superplot(d)
    sc.plot.image(d)
    sc.plot.scatter3d(d)


@pytest.mark.skip(reason="RuntimeError: Only the simple case histograms may "
                  "be constructed for now: 2 dims including events.")
def test_plot_3d_events_data_with_int_bins():
    d = make_events_dataset(ndim=3)
    plot(d, bins=50)


@pytest.mark.skip(reason="RuntimeError: Only the simple case histograms may "
                  "be constructed for now: 2 dims including events.")
def test_plot_3d_events_data_with_nparray_bins():
    d = make_events_dataset(ndim=3)
    plot(d, bins=np.linspace(0.0, 105.0, 50))


@pytest.mark.skip(reason="RuntimeError: Only the simple case histograms may "
                  "be constructed for now: 2 dims including events.")
def test_plot_3d_events_data_with_Variable_bins():
    d = make_events_dataset(ndim=3)
    bins = sc.Variable(['tof'],
                       values=np.linspace(0.0, 105.0, 50),
                       unit=sc.units.us)
    plot(d, bins=bins)


def test_plot_variable_3d():
    N = 50
    v3d = sc.Variable(['tof', 'x', 'y'],
                      values=np.random.rand(N, N, N),
                      unit=sc.units.m)
    plot(v3d)


def test_plot_realigned_3d():
    d = make_events_dataset(ndim=2)
    tbins = sc.Variable(dims=['tof'], unit=sc.units.us, values=np.arange(100.))
    r = sc.realign(d, {'tof': tbins})
    plot(r)


def test_plot_4d_with_masks():
    data = sc.DataArray(data=sc.Variable(
        dims=['pack', 'tube', 'straw', 'pixel'],
        values=np.random.rand(2, 8, 7, 256)),
                        coords={})
    a = np.sin(np.linspace(0, 3.14, num=256))
    data += sc.Variable(dims=['pixel'], values=a)
    data.masks['tube_ends'] = sc.Variable(dims=['pixel'],
                                          values=np.where(
                                              a > 0.5, True, False))
    plot(data)


def test_plot_4d_with_masks_projection_3d():
    data = sc.DataArray(data=sc.Variable(
        dims=['pack', 'tube', 'straw', 'pixel'],
        values=np.random.rand(2, 8, 7, 256)),
                        coords={})
    a = np.sin(np.linspace(0, 3.14, num=256))
    data += sc.Variable(dims=['pixel'], values=a)
    data.masks['tube_ends'] = sc.Variable(dims=['pixel'],
                                          values=np.where(
                                              a > 0.5, True, False))
    plot(data, projection="3d")


def test_plot_3d_data_with_ragged_bins():
    N = 10
    M = 8
    L = 5
    x = np.arange(N + 1).astype(np.float64)
    y = np.arange(M).astype(np.float64)
    z = np.arange(L).astype(np.float64)
    zz, yy, xx = np.meshgrid(z, y, x, indexing='ij')
    a = np.random.random([L, M, N])
    for i in range(M):
        for j in range(L):
            xx[j, i, :] *= (i + j + 1.0)
    d = sc.Dataset()
    d.coords['x'] = sc.Variable(['z', 'y', 'x'], values=xx, unit=sc.units.m)
    d.coords['y'] = sc.Variable(['y'], values=y, unit=sc.units.m)
    d.coords['z'] = sc.Variable(['z'], values=z, unit=sc.units.m)
    d['a'] = sc.Variable(['z', 'y', 'x'], values=a, unit=sc.units.counts)
    plot(d)
