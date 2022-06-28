# Invert a PROJ projection using optimization.
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
from typing import Optional
from .cdll import _cdll
from ctypes import c_double, POINTER, c_char_p, c_ulong


def invert_proj(proj_str: str, x: np.ndarray, y: np.ndarray,
                lon0: Optional[float] = None,
                lat0: Optional[float] = None) -> np.ndarray:
    """
    
    """
    # Sanity:
    x = np.ascontiguousarray(x, dtype=np.double)
    y = np.ascontiguousarray(y, dtype=np.double)
    proj_str = str(proj_str)

    N = x.size
    if y.size != N:
        raise RuntimeError("Size of x and y must be equal.")

    # Determine the projection
    proj_split = [p.split("=") for p in proj_str.split()]
    strip_plus = lambda s : s[1:] if len(s) > 0 and s[0] == '+' else s
    params = {strip_plus(p[0]) : p[1] for p in proj_split if len(p) > 1}
    if "proj" not in params:
        raise RuntimeError("No projection given.")

    # Create the output buffer and call C code:
    lola = np.zeros((N,2))
    proj_str = proj_str.encode("ascii")
    res = _cdll.inverse_project_data_optimize(c_char_p(proj_str), c_ulong(N),
                                          x.ctypes.data_as(POINTER(c_double)),
                                          y.ctypes.data_as(POINTER(c_double)),
                                          lola.ctypes.data_as(POINTER(c_double))
                                        )

    if res != 0:
        raise RuntimeError("Inverse projection failed.")

    lola = lola.T

    return lola[0,:], lola[1,:]

