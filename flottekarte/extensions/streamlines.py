# Compute stream line polygons.
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
from ctypes import c_double, POINTER, c_char_p, c_ulong, c_size_t, byref
from .azimuth import unwrap_azimuth_field


def streamlines(
        xmin: float,
        xmax: float,
        nx: int,
        ymin: float,
        ymax: float,
        ny: int,
        alpha: NDArray[np.float64],
        p0: NDArray[np.float64],
        p1: NDArray[np.float64],
        r: float,
        ds_min: float,
        epsilon: float,
        unwrap_azimuth: bool,
        unwrapping_beta: float,
        unwrapping_Nmax: int
    ) -> list[NDArray[np.float64]]:
    """

    """
    # Sanity:
    xmin_c = c_double(xmin)
    xmax_c = c_double(xmax)
    ymin_c = c_double(ymin)
    ymax_c = c_double(ymax)
    nx_c = c_size_t(nx)
    ny_c = c_size_t(ny)
    shape = (int(nx), int(ny))

    # Smooth the azimuth field by flipping orientations:
    if unwrap_azimuth:
        alpha = unwrap_azimuth_field(
            alpha, unwrapping_beta, unwrapping_Nmax, inplace=False
        )

    # Numpy arrays and underlying buffers. Ensure that they are of the
    # correct shape.
    p0 = np.ascontiguousarray(p0, dtype=np.double)
    p1 = np.ascontiguousarray(p0, dtype=np.double)

    if np.shape(alpha) != shape:
        raise RuntimeError("Shape of alpha must be (nx,ny).")
    if np.shape(p0) != shape:
        raise RuntimeError("Shape of p0 must be (nx,ny).")
    if np.shape(p1) != shape:
        raise RuntimeError("Shape of p1 must be (nx,ny).")

    z = np.empty((*shape, 3))
    z[..., 0] = alpha
    z[..., 1] = p0
    z[..., 2] = p1

    epsilon = float(epsilon)
    if epsilon <= 0.0:
        raise ValueError("epsilon needs to be positive.")

    r_c = c_double(r)
    ds_min_c = c_double(ds_min)
    epsilon_c = c_double(epsilon)

    # Now compute the streamlines:
    struct_id = c_size_t(0)
    res = _cdll.compute_streamlines(
        xmin_c, xmax_c, nx_c, ymin_c, ymax_c, ny_c,
        z.ctypes.data_as(POINTER(c_double)), c_size_t(z.size),
        r_c, ds_min_c, epsilon_c, byref(struct_id)
    )
    if res != 0:
        raise RuntimeError("Computing streamlines failed. Code: "
            + str(int(res))
        )

    # Next copy the polygons from the C++ side:
    Npoly = int(_cdll.get_streamline_polygon_count(struct_id))
    if Npoly == 0:
        print("Npoly == 0")
        return []

    polygons = []
    success = True
    for i in range(Npoly):
        # Figure out how large this polygon is:
        p = c_size_t(i)
        Ni = int(_cdll.get_streamline_polygon_size(struct_id, p))
        if Ni == 0:
            print("N[",i,"] == 0")
            continue

        # The size gives the number of coordinate tuples.
        # Create the correctly shaped buffer here:
        poly_i = np.empty((Ni,2))
        polygons.append(poly_i)

        # Now copy:
        res = int(_cdll.save_streamline_polygon(
            struct_id, p,
            poly_i.ctypes.data_as(POINTER(c_double)),
            c_size_t(poly_i.size)
        ))
        if res != 0:
            print("copy error: res =",res)
            success = False
            break

    # Clear the C++ structures:
    res = _cdll.delete_streamline_struct(struct_id)

    # Now check whether we succeeded:
    if not success:
        raise RuntimeError("Failed to copy some streamline polygons.")

    return polygons