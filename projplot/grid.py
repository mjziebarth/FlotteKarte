# Plotting grid lines.

import numpy as np
from math import sqrt
from .extensions.grid import grid_path
from matplotlib.path import Path
from matplotlib.patches import PathPatch


def plot_grid(ax, xlim: tuple[float,float], ylim: tuple[float,float],
              proj_str: str, tick_spacing_degree: int,
              bisection_offset='auto', minimum_node_distance='auto',
              max_lat=90., linewidth=0.8, **kwargs):
    """
    Plots a geographic coordinate grid.
    """
    diagonal = sqrt((xlim[1] - xlim[0])**2 + (ylim[1] - ylim[0])**2)
    if bisection_offset == 'auto':
        bisection_offset = 1e-2 * diagonal
    if minimum_node_distance == 'auto':
        minimum_node_distance = 1e-2 * diagonal

    # Use backend to generate the path:
    vertices, codes = grid_path(proj_str, *xlim, *ylim, tick_spacing_degree,
                                bisection_offset, minimum_node_distance,
                                max_lat)

    # Matplotlib path and PathPatch:
    path = Path(vertices, codes=codes)
    kwargs["linewidth"] = linewidth
    if "facecolor" not in kwargs:
        kwargs["facecolor"] = 'none'
    ax.add_patch(PathPatch(path, **kwargs))
