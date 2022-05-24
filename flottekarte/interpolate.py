# Routines for interpolating irregulary structured data onto regular grids.
# This file is part of the flottekarte python module.
#
# Authors: Malte J. Ziebarth (ziebarth@gfz-potsdam.de)
#
# Copyright (C) 2022 Deutsches GeoForschungsZentrum Potsdam
#
# Licensed under the EUPL, Version 1.2 or – as soon they will be approved by
# the European Commission - subsequent versions of the EUPL (the "Licence");
# You may not use this work except in compliance with the Licence.
# You may obtain a copy of the Licence at:
#
# https://joinup.ec.europa.eu/collection/eupl/eupl-text-eupl-12
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the Licence is distributed on an "AS IS" basis,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the Licence for the specific language governing permissions and
# limitations under the Licence.

import numpy as np
from math import ceil
from typing import Tuple
from scipy.spatial import KDTree
from scipy.sparse import coo_array
from scipy.interpolate import interp2d, LinearNDInterpolator


def interpolate_to_raster(xyz: np.ndarray, xlim: Tuple[float,float],
                          ylim: Tuple[float,float],
                          resolution: float) -> Tuple[np.ndarray, np.ndarray,
                                                      np.ndarray]:
    """
    Interpolate unstructured data `xyz` to a freshly created raster,
    the latter spanning the `xlim` and `ylim` coordinate ranges with
    spacing between grid points being approximately `resolution`.
    The data is interpolated using linear interpolation (SciPy
    LinearNDInterpolator).

    Returns:
       xg, yg, zg

    Where `xg` and `yg`, shapes (nx,) and (ny,), denote the coordinates
    of the new raster. `zg`, shape (ny,nx), is the interpolated data
    raster.
    """
    interp = LinearNDInterpolator(xyz[:,:2], xyz[:,2])

    # Generate the grid:
    nx = ceil((xlim[1]-xlim[0])/resolution)
    ny = ceil((ylim[1]-ylim[0])/resolution)
    xg = np.linspace(xlim[0], xlim[1], nx)
    yg = np.linspace(ylim[0], ylim[1], ny)

    return xg, yg, interp(xg[np.newaxis,:], yg[:,np.newaxis])


def interpolate_to_raster_gauss(xyz: np.ndarray, xlim: Tuple[float,float],
                                ylim: Tuple[float,float], resolution: float,
                                rscale: float) -> Tuple[np.ndarray, np.ndarray,
                                                        np.ndarray]:
    """
    Interpolate unstructured data `xyz` to a freshly created raster,
    the latter spanning the `xlim` and `ylim` coordinate ranges with
    spacing between grid points being approximately `resolution`.
    The data is interpolated using a Gaussian kernel of bandwidth `rscale`,
    taking into account data points within 3*rscale distance
    (or a minimum of 5).

    Returns:
       xg, yg, zg

    Where `xg` and `yg`, shapes (nx,) and (ny,), denote the coordinates
    of the new raster. `zg`, shape (ny,nx), is the interpolated data
    raster.
    """
    # Generate the grid:
    nx = ceil((xlim[1]-xlim[0])/resolution)
    ny = ceil((ylim[1]-ylim[0])/resolution)
    xg = np.linspace(xlim[0], xlim[1], nx)
    yg = np.linspace(ylim[0], ylim[1], ny)

    xgg,ygg = np.meshgrid(xg,yg)

    # Create query trees:
    tree_src = KDTree(xyz[:,:2])
    xy_dest = np.stack((xgg.flat, ygg.flat), axis=1)
    tree_dest = KDTree(xy_dest)

    ## Interpolate:
    res = tree_dest.query_ball_tree(tree_src, 3*rscale)
    few_neighbors = [len(nb) < 5 for nb in res]
    if any(few_neighbors):
        few_neighbors = np.argwhere(few_neighbors)
        nn = tree_src.query(xy_dest[few_neighbors], 5)
        for j,i in enumerate(few_neighbors):
            res[i] = nn[j]

    rs2 = rscale ** 2
    weights = [np.exp(-0.5*np.sum((xyz[r,:2] - xy_dest[i,:][np.newaxis,:])**2,
                                  axis=1) / rs2)
               for i,r in enumerate(res)]
    weights = [w / w.sum() for w in weights]
    zg = np.zeros((ny,nx))
    zg.flat[:] = [(w * xyz[r,2]).sum() for w,r in zip(weights,res)]

    return xg, yg, zg 
