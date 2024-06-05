# A two-dimensional tensor field representation
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
from pyproj import Proj
from numpy.typing import NDArray
from .extensions.azimuth import azimuth_geographic_to_local_on_grid,\
    unwrap_azimuth_field

class TensorField2D:
    """
    A two-dimensional tensor field class.
    """

    xmin: float
    xmax: float
    nx: int
    ymin: float
    ymax: float
    ny: int
    # A note on the angle definition we use here.
    # We define 'alpha' to be the clockwise angle from the y axis
    # with x pointing eastward. I.e. when using a Plate Carée, 0°
    # points north and 90° points east.
    # Also, we save the angle in radians.
    alpha_rad: NDArray
    p0: NDArray
    p1: NDArray

    unwrapped: bool

    def __init__(self,
            xmin: float,
            xmax: float,
            nx: int,
            ymin: float,
            ymax: float,
            ny: int,
            p0: NDArray,
            p1: NDArray,
            *args,
            geographic_azimuth: NDArray[np.float64] | None = None,
            projected_azimuth: NDArray[np.float64] | None = None,
            proj_str: str | None = None,
            stencil_delta: float = 1e-5,
            unwrap_azimuth: bool = False,
            unwrapping_beta: float = 1.5
        ):
        self.xmin = float(xmin)
        self.xmax = float(xmax)
        self.nx = int(nx)
        self.ymin = float(ymin)
        self.ymax = float(ymax)
        self.ny = int(ny)
        shape = (nx, ny)
        if geographic_azimuth is None:
            if projected_azimuth is None:
                raise ValueError("One of 'geographic_azimuth' or "
                    "'projected_azimuth' needs to be provided."
                )
            alpha = np.deg2rad(projected_azimuth)
        else:
            if projected_azimuth is not None:
                raise ValueError("Only one of 'geographic_azimuth' or "
                    "'projected_azimuth' can be provided."
                )
            if proj_str is None:
                raise ValueError("If geographic_azimuth is given, "
                    "the proj_str has to be given to allow conversion to "
                    "projected azimuth."
                )

            # Project the geographic azimuth:
            alpha = azimuth_geographic_to_local_on_grid(
                proj_str, xmin, xmax, nx, ymin, ymax, ny,
                np.deg2rad(geographic_azimuth), inplace=True,
                stencil_delta = float(stencil_delta)
            )

        # If wished, unwrap the azimuth field:
        if unwrap_azimuth:
            alpha = unwrap_azimuth_field(
                alpha, unwrapping_beta, inplace=True
            )
            self.unwrapped = True
        else:
            self.unwrapped = False


        self.alpha_rad = alpha
        self.p0 = np.ascontiguousarray(p0, dtype=np.double)
        self.p1 = np.ascontiguousarray(p1, dtype=np.double)
        if self.alpha_rad.shape != shape or self.p0.shape != shape \
                or self.p1.shape != shape:
            raise ValueError("All of 'p0', 'p1', and the azimuth need to be "
                "of shape (" + str(nx) + "," + str(ny) + ")"
            )
