# The Map class.
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
from typing import Union
from pyproj import Proj


# Import of functionality from other source files of this package:
from .axes import generate_axes_grid
from .grid import plot_grid
from .data import GeoJSON
from .extensions.boundary import map_boundary
from .extent import automatic_map_extents


class Map:
    """
    A minimalistic map plotting tool using PROJ, pyproj,
    and Matplotlib.
    """
    def __init__(self, proj_str, ax, xlim=None, ylim=None, proj=None,
                 bisection_offset='auto', minimum_node_distance='auto',
                 atol='auto'):
        self.proj_str = proj_str
        self.ax = ax
        if proj is None:
            proj = Proj(proj_str)

        # Map limits, potentially automatically generated:
        if xlim is None or ylim is None:
            xmin_a,xmax_a,ymin_a,ymax_a = automatic_map_extents(proj_str, proj)
            if xlim is None:
                xlim = xmin_a, xmax_a
            if ylim is None:
                ylim = ymin_a, ymax_a
        self.xlim = xlim
        self.ylim = ylim

        self.proj = proj
        self._grid_cuts = []
        self._bisection_offset = bisection_offset
        self._minimum_node_distance = minimum_node_distance
        self._atol = atol

        # Compute the map boundary:
        self.boundary_vertices, self.boundary_normal_angles \
            = map_boundary(proj_str, *xlim, *ylim, atol, bisection_offset,
                           minimum_node_distance)

        # Some initialization:
        ax.set_aspect('equal')
        ax.set_axis_off()

    @staticmethod
    def for_data(lon: np.ndarray, lat:np.ndarray, proj_str: str,
                 ax, margin_rel: float = 0.05):
        """
        Generate a Map fitting some data in a given projection.
        """
        proj = Proj(proj_str)
        x,y = proj(lon, lat)
        xlim = (x.min(), x.max())
        ylim = (y.min(), y.max())
        Dx = xlim[1] - xlim[0]
        Dy = ylim[1] - ylim[0]
        margin = margin_rel * max(Dx,Dy)
        return Map(proj_str, ax, (xlim[0]-margin, xlim[1]+margin),
                   (ylim[0]-margin, ylim[1]+margin), proj)

    def plot_axes(self, tick_spacing: float = 1.0,
                  which_ticks: Union[str,list] = 'auto',
                  linewidth: float = 0.8, fontsize: float = 8.0,
                  rotate_labels: bool = True, ticks_bot: bool = True,
                  ticks_top: bool = True, ticks_left: bool = True,
                  ticks_right: bool = True,
                  remove_overlaps: bool = True):
        """
        Add axes and ticks to the plot.
        """
        generate_axes_grid(self.ax, self.boundary_vertices,
                           self.boundary_normal_angles, self.proj_str,
                           tick_spacing=tick_spacing, which_ticks=which_ticks,
                           proj=self.proj, fontsize=fontsize,
                           grid_cuts=self._grid_cuts,
                           rotate_labels=rotate_labels, ticks_bot=ticks_bot,
                           ticks_top=ticks_top, ticks_left=ticks_left,
                           ticks_right=ticks_right,
                           remove_overlaps=remove_overlaps)

    def plot_grid(self, tick_spacing_degree: int = 10, max_lat: float = 90.0,
                  bisection_offset=None, minimum_node_distance=None,
                  linewidth=0.8, **kwargs):
        """
        Add a geographic coordinte grid to the plot.
        """
        if bisection_offset is None:
            bisection_offset = self._bisection_offset
        if minimum_node_distance is None:
            minimum_node_distance = self._minimum_node_distance

        cuts, cut_tick_type, cut_coordinates = \
           plot_grid(self.ax, self.xlim, self.ylim, self.proj_str,
                     tick_spacing_degree, bisection_offset=bisection_offset,
                     minimum_node_distance=minimum_node_distance,
                     max_lat=max_lat, linewidth=linewidth, **kwargs)
        self._grid_cuts = [(float(x), float(y), "lon" if tt == 0 else "lat",
                            float(coordinate))
                           for (x,y), tt, coordinate in zip(cuts, cut_tick_type,
                                                            cut_coordinates)]

    def add_data(self, dataset, **kwargs):
        """
        Add data loaded from the `data` submodule.
        """
        if isinstance(dataset,GeoJSON):
            polygons = dataset.get_polygon_patches(**kwargs)
            self.ax.add_collection(polygons)
        else:
            raise TypeError("`dataset` must be a GeoJSON object.")
