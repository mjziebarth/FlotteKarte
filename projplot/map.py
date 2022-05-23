# The Map class.

import numpy as np
from pyproj import Proj


# Import of functionality from other source files of this package:
from .axes import generate_axes_grid
from .grid import plot_grid


class Map:
    """
    A minimalistic map plotting tool using PROJ, pyproj,
    and Matplotlib.
    """
    def __init__(self, proj_str, ax, xlim, ylim, proj=None):
        self.proj_str = proj_str
        self.ax = ax
        self.xlim = xlim
        self.ylim = ylim
        if proj is None:
            proj = Proj(proj_str)
        self.proj = proj

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
                  tick_bot: str = 'lon', tick_top: str = 'lon',
                  tick_left:str = 'lat', tick_right: str = 'lat',
                  linewidth: float = 0.8, fontsize: float = 8.0):
        """
        Add axes and ticks to the plot.
        """
        generate_axes_grid(self.ax, self.xlim, self.ylim, self.proj_str,
                           tick_spacing=tick_spacing, tick_bot=tick_bot,
                           tick_top=tick_top, tick_left=tick_left,
                           tick_right=tick_right, proj=self.proj,
                           fontsize=fontsize)

    def plot_grid(self, tick_spacing_degree: int = 10, max_lat: float = 90.0,
                  bisection_offset='auto', minimum_node_distance='auto',
                  cut_angle_at_degrees: float = 90.0,
                  linewidth=0.8, **kwargs):
        """
        Add a geographic coordinte grid to the plot.
        """
        plot_grid(self.ax, self.xlim, self.ylim, self.proj_str,
                  tick_spacing_degree, bisection_offset=bisection_offset,
                  minimum_node_distance=minimum_node_distance,
                  max_lat=max_lat, cut_angle_at_degrees=cut_angle_at_degrees,
                  linewidth=linewidth, **kwargs)
