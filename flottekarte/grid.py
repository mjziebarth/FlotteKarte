# Plotting grid lines.
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
from typing import Tuple, Union, Literal
from math import sqrt
from .extensions.grid import grid_path
from matplotlib.path import Path
from matplotlib.patches import PathPatch


def plot_grid(ax, xlim: Tuple[float,float], ylim: Tuple[float,float],
              proj_str: str, tick_spacing_degree: int,
              bisection_offset: Union[Literal['auto'],float] = 'auto',
              minimum_node_distance: Union[Literal['auto'],float] = 'auto',
              max_lat=90., linewidth=0.8, **kwargs):
    """
    Plots a geographic coordinate grid.
    """
    diagonal = sqrt((xlim[1] - xlim[0])**2 + (ylim[1] - ylim[0])**2)
    if bisection_offset == 'auto':
        bisection_offset = 1e-2 * diagonal
    if minimum_node_distance == 'auto':
        minimum_node_distance = 1e-3 * diagonal

    # Use backend to generate the path:
    vertices, codes, cuts, cut_tick_type, cut_coordinates \
       = grid_path(proj_str, *xlim, *ylim, tick_spacing_degree,
                   bisection_offset, minimum_node_distance, max_lat)

    # Matplotlib path and PathPatch:
    path = Path(vertices, codes=codes)
    kwargs["linewidth"] = linewidth
    if "facecolor" not in kwargs:
        kwargs["facecolor"] = 'none'
    ax.add_patch(PathPatch(path, **kwargs))

    return cuts, cut_tick_type, cut_coordinates
