# Compute a path of  of coordinates.
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
from .cdll import _cdll
from typing import Tuple
from ctypes import c_double, c_int, POINTER, c_char_p, c_ulong, c_void_p, \
                   byref, c_ubyte


def grid_path(proj_str: str, xmin: float, xmax: float, ymin: float, ymax: float,
              tick_spacing_degree: int, bisection_offset: float,
              minimum_node_distance: float,
              max_lat: float) -> Tuple[np.ndarray,np.ndarray,np.ndarray]:
    """
    Computes a Matplotlib path comprising the lines of a geograpnic coordinate
    grid.

    Returns:
       vertices, codes
    """
    # Sanity:
    proj_str = str(proj_str)
    xmin = float(xmin)
    xmax = float(xmax)
    if xmin >= xmax:
        raise RuntimeError("`xmin` has to be smaller than `xmax`.")
    ymin = float(ymin)
    ymax = float(ymax)
    if ymin >= ymax:
        raise RuntimeError("`ymin` has to be smaller than `ymax`.")
    tick_spacing_degree = int(tick_spacing_degree)
    if tick_spacing_degree == 0:
        raise RuntimeError("`tick_spacing_degree` may not be zero.")
    bisection_offset = float(bisection_offset)
    if bisection_offset <= 0:
        raise RuntimeError("`bisection_offset` has to be positive.")
    minimum_node_distance = float(minimum_node_distance)
    if minimum_node_distance <= 0:
        raise RuntimeError("`minimum_node_distance` has to be positive.")
    max_lat = float(max_lat)
    if max_lat <= 0 or max_lat > 90.0:
        raise RuntimeError("`max_lat` has to be positive and may not be larger "
                           "than 90.")

    # Determine the projection
    proj_split = [p.split("=") for p in proj_str.split()]
    strip_plus = lambda s : s[1:] if len(s) > 0 and s[0] == '+' else s
    params = {strip_plus(p[0]) : p[1] for p in proj_split if len(p) > 1}
    if "proj" not in params:
        raise RuntimeError("No projection given.")

    # Create the output buffers and call C code:
    struct_ptr = c_void_p(0)
    Npath = c_ulong(0)
    Ncut = c_ulong(0)
    proj_str = proj_str.encode("ascii")

    res = _cdll.compute_grid_lines(c_char_p(proj_str), c_double(xmin),
                                   c_double(xmax), c_double(ymin),
                                   c_double(ymax), c_int(tick_spacing_degree),
                                   c_double(bisection_offset),
                                   c_double(minimum_node_distance),
                                   c_double(max_lat),
                                   byref(struct_ptr), byref(Npath),
                                   byref(Ncut));

    if res != 0:
        # Clean up:
        print("struct_ptr:", struct_ptr)
        if struct_ptr:
            _cdll.clean_grid_lines_struct(struct_ptr)
        raise RuntimeError("compute_grid_lines backend failed. Error code "
                           + str(res))

    # Create the path array:
    vertices = np.zeros((Npath.value, 2))
    codes = np.zeros(Npath.value, dtype=np.uint8)

    # Cut array:
    cuts = np.zeros((Ncut.value, 2))
    cut_tick_type = np.zeros(Ncut.value, dtype=np.uint8)
    cut_coordinates = np.zeros(Ncut.value)

    # Fill the path to the numpy arrays:
    res = _cdll.save_grid_lines(struct_ptr,
                                vertices.ctypes.data_as(POINTER(c_double)),
                                codes.ctypes.data_as(POINTER(c_ubyte)),
                                cuts.ctypes.data_as(POINTER(c_double)),
                                cut_tick_type.ctypes.data_as(POINTER(c_double)),
                                cut_coordinates.ctypes
                                               .data_as(POINTER(c_double)))

    # Clean the C++ structures:
    _cdll.clean_grid_lines_struct(struct_ptr)

    if res != 0:
        raise RuntimeError("save_grid_lines backend failed. Error code "
                           + str(res))

    return vertices, codes, cuts, cut_tick_type, cut_coordinates
