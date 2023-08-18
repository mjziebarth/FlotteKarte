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
from typing import Union, Tuple, Optional, Literal
from pyproj import Proj
from matplotlib.axes import Axes


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

    Parameters
    ----------
    proj_str : str
       The legacy PROJ string definition of a projection.
    ax : matplotlib.axes.Axes
       The axes to transform into a map.
    xlim : tuple[float,float], optional
       The maximum *x* extent of the map (:math:`x_\\mathrm{min}`,
       :math:`x_\\mathrm{max}`).
    ylim : tuple[float,float], optional
       The maximum *y* extent of the map (:math:`y_\\mathrm{min}`,
       :math:`y_\\mathrm{max}`).
    proj : pyproj.Proj, optional
       A :class:`pyproj.Proj` instance that corresponds to `proj_str`.

    bisection_offset : float, optional
       Critical offset that determines whether a path segment will
       be refined. If the center of a spherical geodesic is offset
       further than this value from the center of the line segment
       in projected coordinates, the line segment is split (bisected).
       Given in projected coordinates. If not given, will be set to
       a fraction of the map extent.
    minimum_node_distance : float, optional
       Stopping criterion for the segment refinement. If the two
       end points of a line segment are closer than this distance,
       the segment will not be split in refinement.
       Given in projected coordinates. If not given, will be set
       to a fraction of the map extent.
    atol : float, optional
       Absolute tolerance to which a non-rectangular map boundary will
       be determined. This parameter is used when it has been
       determined that the bounding rectangle defined by `xlim` and
       `ylim` contains space that is not actually part of the image
       of the projection (e.g. when working with a global projection
       or close to the boundaries of an oblique projection).
       Given in projected coordinates. If not given, will be set to
       a fraction of the map extent.
    """
    def __init__(self, proj_str: str, ax: Axes,
                 xlim: Optional[Tuple[float,float]] = None,
                 ylim: Optional[Tuple[float,float]] = None,
                 proj: Optional[Proj] = None,
                 bisection_offset: Union[Literal['auto'],float] ='auto',
                 minimum_node_distance: Union[Literal['auto'],float] = 'auto',
                 atol: Union[Literal['auto'],float] = 'auto'):
        if 'type=crs' in proj_str:
            raise ValueError("Likely running into PROJ error if 'type=crs' "
                             "is given. Please remove this argument.")
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
    def for_data(lon: np.ndarray, lat: np.ndarray, proj_str: str,
                 ax, margin_rel: float = 0.05):
        """
        Generate a Map fitting some data in a given projection.

        Parameters
        ----------
        lon : array_like
           Longitudes of the data set.
        lat : array_like
           Latitudes of the data set.
        proj_str : str
           The legacy PROJ string definition of a projection.
        margin_rel : float, optional
           Introduce a margin relative to the extent of the
           data in projected coordinates.


        Returns
        -------
        :class:`Map`
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

        Parameters
        ----------
        tick_spacing : float, optional
           Spacing of ticks in degrees.
        which_ticks : str or list, optional
           The type of ticks to favor, longitude or latitude.
           If given as a ``list``, has to be a four-element
           list that contains one of ``lon`` or ``lat`` per
           part of the boundary.
           **Note: currently not used.**
        linewidth : float, optional
           Line width of boundary and ticks.
        fontsize : float, optional
           Font size, in point, of coordinate tick labels.
        rotate_labels : bool, optional
           Whether to rotate the labels parallel to the
           map boundary at the tick position.
        ticks_bot : bool, optional
           Whether to plot ticks on the bottom axis. For a
           non-rectangular map shape, 'bottom axis' refers
           to all ticks that extend beyond ``ymin``.
        ticks_top : bool, optional
           Whether to plot ticks on the top axis. For a
           non-rectangular map shape, 'top axis' refers
           to all ticks that extend beyond ``ymax``.
        ticks_left : bool, optional
           Whether to plot ticks on the left axis. For a
           non-rectangular map shape, 'left axis' refers
           to all ticks that extend beyond ``xmin``.
        ticks_right : bool, optional
           Whether to plot ticks on the right axis. For a
           non-rectangular map shape, 'right axis' refers
           to all ticks that extend beyond ``xmax``.
        remove_overlaps : bool, optional
           Whether to apply the overlap removal algorithm.


        Notes
        -----
        The overlap algorithm uses a tick weighting scheme
        to remove the least weighted ticks whenever two ticks
        overlap. This tick weighting scheme considers, in
        descending order of importance, the following
        attributes of each tick:

        1) Is the tick requested?
        2) Is the tick part of a grid line intersecting the
           boundary?
        3) Is the tick at a "nice" value? (E.g. 0°, 90°,
           zero modul for a number of common divisiors)

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
                  bisection_offset: Optional[Union[str,float]] = None,
                  minimum_node_distance: Optional[Union[str,float]] = None,
                  linewidth: float = 0.8, **kwargs):
        """
        Add a geographic coordinte grid to the plot.

        Parameters
        ----------
        tick_spacing_degree : int, optional
           Spacing of ticks in degrees.
        max_lat : float, optional
           Maximum absolute latitude up to which to draw the
           grid lines.
        bisection_offset : float, optional
           Critical offset that determines whether a path segment will
           be refined. If the center of a spherical geodesic is offset
           further than this value from the center of the line segment
           in projected coordinates, the line segment is split (bisected).
           Given in projected coordinates. If not given, will be set to
           a fraction of the map extent.
        minimum_node_distance : float, optional
           Stopping criterion for the segment refinement. If the two
           end points of a line segment are closer than this distance,
           the segment will not be split in refinement.
           Given in projected coordinates. If not given, will be set
           to a fraction of the map extent.
        linewidth : float, optional
           Width of the grid lines.
        kwargs : optional
           Passed to :py:class:`matplotlib.LineCollection`.
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

    def add_data(self, dataset: Union[GeoJSON], **kwargs):
        """
        Add data loaded from the ``data`` submodule.

        Parameters
        ----------
        dataset : GeoJSON
           A data set from the `flottekarte.data` submodule.
           Currently that includes only the :class:`GeoJSON` class.

        See Also
        --------
        GeoJSON
        """
        if isinstance(dataset,GeoJSON):
            polygons = dataset.get_polygon_patches(**kwargs)
            self.ax.add_collection(polygons)
        else:
            raise TypeError("`dataset` must be a GeoJSON object.")
