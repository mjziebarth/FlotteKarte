# Axes ticks computation.
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
from ctypes import c_double, POINTER, c_char_p, c_uint, c_ubyte, c_ulong, byref
from typing import Optional, Tuple

def compute_axes_ticks(proj_str: str, boundary: np.ndarray,
                       tick_spacing: float) \
    -> Tuple[np.ndarray,np.ndarray,np.ndarray]:
    """

    """
    # Sanity:
    boundary = np.ascontiguousarray(boundary)
    if boundary.ndim != 2 or  boundary.shape[1] != 2:
        raise RuntimeError("`boundary` must have shape (Nx2)")
    tick_spacing = float(tick_spacing)
    proj_str = str(proj_str)

    # Determine the projection
    proj_split = [p.split("=") for p in proj_str.split()]
    strip_plus = lambda s : s[1:] if len(s) > 0 and s[0] == '+' else s
    params = {strip_plus(p[0]) : p[1] for p in proj_split if len(p) > 1}
    if "proj" not in params:
        raise RuntimeError("No projection given.")

    # Create the output buffers and call C code:
    Nmax = 1000
    tick_vertices = np.zeros((Nmax,2))
    segments = np.zeros(Nmax, dtype=np.uintc)
    which_ticks = np.zeros(Nmax, dtype=np.uint8)
    proj_str = proj_str.encode("ascii")
    Nticks = c_uint(0)

    res = _cdll.compute_axes_ticks(c_char_p(proj_str),
                                c_ulong(boundary.shape[0]),
                                boundary.ctypes.data_as(POINTER(c_double)),
                                c_double(tick_spacing), c_uint(Nmax),
                                segments.ctypes.data_as(POINTER(c_uint)),
                                tick_vertices.ctypes.data_as(POINTER(c_double)),
                                which_ticks.ctypes.data_as(POINTER(c_ubyte)),
                                byref(Nticks))

    Nticks = Nticks.value

    if res != 0:
        raise RuntimeError("grad_east_north backend failed.")

    return (tick_vertices[:Nticks].copy(), which_ticks[:Nticks].copy(),
            segments[:Nticks].copy())
