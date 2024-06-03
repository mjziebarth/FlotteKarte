# Conversion of vector orientations (azimuth) between the physical/geographic
# space and the projected space.
#
# Authors: Malte J. Ziebarth (malte.ziebarth@tum.de)
#
# Copyright (C) 2024 Technische Universität München
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
from numpy.typing import NDArray
from typing import Optional
from .cdll import _cdll
from ctypes import c_double, c_size_t, c_char_p, POINTER


def azimuth_geographic_to_local_on_grid(
        proj_str: str,
        xmin: float, xmax: float, nx: int,
        ymin: float, ymax: float, ny: int,
        azimuth_rad: NDArray[np.float64],
        inplace: bool = False,
        stencil_delta: float = 1e-5
    ) -> NDArray[np.float64]:
    """
    Computes the local azimuth from the geographic azimuth.
    """
    proj_str = str(proj_str).encode("ascii")
    proj_str_c = c_char_p(proj_str)
    xmin_c = c_double(xmin)
    xmax_c = c_double(xmax)
    ymin_c = c_double(ymin)
    ymax_c = c_double(ymax)
    nx_c = c_size_t(nx)
    ny_c = c_size_t(ny)
    shape = (int(nx), int(ny))

    # Numpy arrays and underlying buffers. Ensure that they are of the
    # correct shape.
    alpha = np.array(
        azimuth_rad, dtype=np.double,
        order='C', copy = not inplace
    )

    if alpha.shape != shape:
        raise RuntimeError("Shape of azimuth must be (nx,ny).")


    # Now compute the azimuth inplace:
    struct_id = c_size_t(0)
    res = _cdll.azimuth_geographic_to_local_on_grid_inplace(
        proj_str_c,
        xmin_c, xmax_c, nx_c, ymin_c, ymax_c, ny_c,
        alpha.ctypes.data_as(POINTER(c_double)), c_size_t(alpha.size),
        c_double(stencil_delta)
    )

    if res != 0:
        raise RuntimeError("Did not succeed to compute local azimuth.")

    return alpha