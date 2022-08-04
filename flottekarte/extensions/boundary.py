# Compute the map boundary.
#
# Authors: Malte J. Ziebarth (ziebarth@gfz-potsdam.de)
#
# Copyright (C) 2022 Deutsches GeoForschungsZentrum Potsdam,
#                    Malte J. Ziebarth
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
from math import sqrt
from .cdll import _cdll
from typing import Tuple, Union
from ctypes import c_double, c_int, POINTER, c_char_p, c_ulong, c_void_p, \
                   byref, c_ubyte


def map_boundary(proj_str: str, xmin: float, xmax: float, ymin: float,
                 ymax: float, atol: Union[float,str],
                 bisection_offset: Union[float,str],
                 minimum_node_distance: Union[float,str]) \
    -> Tuple[np.ndarray,np.ndarray]:
    """
    Boundary of a map from a projection and a bounding box.

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

    diagonal = sqrt((xmax - xmin)**2 + (ymax - ymin)**2)
    if bisection_offset == 'auto':
        bisection_offset = 1e-3 * diagonal
    else:
        bisection_offset = float(bisection_offset)
    if minimum_node_distance == 'auto':
        minimum_node_distance = 1e-3 * diagonal
    else:
        minimum_node_distance = float(minimum_node_distance)
    if atol == 'auto':
        atol = 1e-4 * diagonal
    else:
        atol = float(atol)

    # Determine the projection
    proj_split = [p.split("=") for p in proj_str.split()]
    strip_plus = lambda s : s[1:] if len(s) > 0 and s[0] == '+' else s
    params = {strip_plus(p[0]) : p[1] for p in proj_split if len(p) > 1}
    if "proj" not in params:
        raise RuntimeError("No projection given.")

    # Create the output buffers and call C code:
    struct_ptr = c_void_p(0)
    Nvert = c_ulong(0)
    proj_str = proj_str.encode("ascii")

    res = _cdll.compute_bounding_polygon(c_char_p(proj_str), c_double(xmin),
                                   c_double(xmax), c_double(ymin),
                                   c_double(ymax), c_double(atol),
                                   c_double(bisection_offset),
                                   c_double(minimum_node_distance),
                                   byref(struct_ptr), byref(Nvert));


    if res != 0:
        # Clean up:
        print("struct_ptr:", struct_ptr)
        if struct_ptr:
            _cdll.clean_bounding_polygon_struct(struct_ptr)
        raise RuntimeError("compute_bounding_polygon backend failed. Error "
                           "code " + str(res))

    # Create the polygon and angle arrays:
    vertices = np.zeros((Nvert.value, 2))
    angles = np.zeros(Nvert.value)

    # Fill the path to the numpy arrays:
    res = _cdll.save_bounding_polygon(struct_ptr,
                                vertices.ctypes.data_as(POINTER(c_double)),
                                angles.ctypes.data_as(POINTER(c_double))
    )

    # Clean the C++ structures:
    _cdll.clean_bounding_polygon_struct(struct_ptr)

    if res != 0:
        raise RuntimeError("save_bounding_polygon backend failed. Error code "
                           + str(res))

    return vertices, angles
