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
from typing import Literal
from .cdll import _cdll
from ctypes import c_double, c_size_t, c_char_p, c_uint32, POINTER


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


def unwrap_azimuth_field(
        azimuth_rad: NDArray[np.float64],
        beta: float = 1.5,
        Nmax: int | Literal['auto'] = 'auto',
        inplace: bool = False
    ):
    """
    Attempt to improve the quality of the azimuth field in terms
    of a continuous vector field. This function aims to tackle the
    issues appearing when defining azimuth only over a 180° range
    instead of the full 360° range, and thereafter integrating over
    the resulting vector field.
    """

    # Underlying working buffer. Ensure that it is of the correct shape
    # and continuity:
    alpha = np.array(
        azimuth_rad, dtype=np.double,
        order='C', copy = not inplace
    )
    if alpha.ndim != 2:
        raise TypeError("alpha has to be two-dimensional.")

    if Nmax == 'auto':
        Nmax = alpha.size

    nx_c = c_uint32(alpha.shape[0])
    ny_c = c_uint32(alpha.shape[1])
    Nmax_c = c_size_t(Nmax)
    beta_c = c_double(beta)

    res = _cdll.unwrap_azimuth_field(
        alpha.ctypes.data_as(POINTER(c_double)), nx_c, ny_c, Nmax_c,
        beta_c
    )
    if res != 0:
        raise RuntimeError("The C++ backend of "
            "unwrap_azimuth_field encountered a runtime error."
        )

    return azimuth_rad