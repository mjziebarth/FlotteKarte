# Invert a PROJ projection using optimization.

import pathlib
import numpy as np
from typing import Optional
from ctypes import CDLL, c_double, c_size_t, c_uint, POINTER, c_ushort,\
                   c_char_p, c_ulong

# Load the shared library:
_cdll_path = pathlib.Path(__file__).parent / 'libinverse.so'
_cdll = CDLL(_cdll_path)


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
    projection = params["proj"]

    # Try to estimate some defaults:
    if lon0 is None:
        if "lonc" in params:
            lon0 = float(params["lonc"])
        elif "lon_0" in params:
            lon0 = float(params["lon_0"])
        elif "lon_1" in params:
            lon0 = float(params["lon_1"])
        else:
            lon0 = 0.0
    else:
        lon0 = float(lon0)
    if lat0 is None:
        if "lat_0" in params:
            lat0 = float(params["lat_0"])
        elif "lat_1" in params:
            lat0 = float(params["lat_1"])
        elif "lat_ts" in params:
            lat0 = float(params["lat_ts"])
        else:
            lat0 = 0.0
    else:
        lat0 = float(lat0)

    # Create the output buffer and call C code:
    lola = np.zeros((N,2))
    proj_str = proj_str.encode("ascii")
    _cdll.inverse_project_data_optimize(c_char_p(proj_str), c_ulong(N),
                                        x.ctypes.data_as(POINTER(c_double)),
                                        y.ctypes.data_as(POINTER(c_double)),
                                        c_double(lon0), c_double(lat0),
                                        lola.ctypes.data_as(POINTER(c_double))
                                      )

    lola = lola.T

    return lola[0,:], lola[1,:]

