# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2019 Scipp contributors (https://github.com/scipp)
# @author Simon Heybrock
import numpy as np
import scipp as sc
from scipp import Dim


def make_dataarray(dim1=Dim.X, dim2=Dim.Y, seed=None):
    if seed is not None:
        np.random.seed(seed)
    return sc.DataArray(
        data=sc.Variable(dims=[dim1, dim2], values=np.random.rand(2, 3)),
        coords={
            dim1: sc.Variable([dim1], values=np.arange(2.0), unit=sc.units.m),
            dim2: sc.Variable([dim2], values=np.arange(3.0), unit=sc.units.m)
        },
        labels={'aux': sc.Variable([dim2], values=np.random.rand(3))},
        attrs={'meta': sc.Variable([dim2], values=np.arange(3))})


def test_init():
    d = sc.DataArray(
        data=sc.Variable(dims=[sc.Dim.X], values=np.arange(3)),
        coords={
            sc.Dim.X:
            sc.Variable([sc.Dim.X], values=np.arange(3), unit=sc.units.m),
        },
        labels={'lib1': sc.Variable([sc.Dim.X], values=np.random.rand(3))},
        attrs={'met1': sc.Variable([sc.Dim.X], values=np.arange(3))},
        masks={
            'mask1': sc.Variable([sc.Dim.X], values=np.ones(3, dtype=np.bool))
        })
    assert len(d.attrs) == 1
    assert len(d.labels) == 1
    assert len(d.masks) == 1


def test_in_place_binary_with_variable():
    a = sc.DataArray(
        data=sc.Variable([Dim.X], values=np.arange(10.0)),
        coords={Dim.X: sc.Variable([Dim.X], values=np.arange(10.0))})
    copy = a.copy()

    a += 2.0 * sc.units.dimensionless
    a *= 2.0 * sc.units.m
    a -= 4.0 * sc.units.m
    a /= 2.0 * sc.units.m
    assert a == copy


def test_in_place_binary_with_scalar():
    a = sc.DataArray(data=sc.Variable([Dim.X], values=[10]),
                     coords={Dim.X: sc.Variable([Dim.X], values=[10])})
    copy = a.copy()

    a += 2
    a *= 2
    a -= 4
    a /= 2
    assert a == copy


def test_rename_dims():
    d = make_dataarray(Dim.X, Dim.Y, seed=0)
    d.rename_dims({Dim.Y: Dim.Z})
    assert d == make_dataarray(Dim.X, Dim.Z, seed=0)
    d.rename_dims(dims_dict={Dim.X: Dim.Y, Dim.Z: Dim.X})
    assert d == make_dataarray(Dim.Y, Dim.X, seed=0)


def test_setitem_works_for_view_and_array():
    a = make_dataarray(Dim.X, Dim.Y, seed=0)
    a[Dim.X, :][Dim.X, 0] = a[Dim.X, 1]
    a[Dim.X, 0] = a[Dim.X, 1]
