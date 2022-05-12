# Compute gradients of coordinates.

import numpy as np
from .cdll import _cdll
from ctypes import c_double, POINTER, c_char_p, c_ulong

def gradients_east_north(proj_str: str, lon: np.ndarray, lat: np.ndarray,
                         delta: float = 1e-5) -> tuple[np.ndarray,np.ndarray]:
    """

    """
    # Sanity:
    lon = np.ascontiguousarray(lon, dtype=np.double)
    lat = np.ascontiguousarray(lat, dtype=np.double)
    proj_str = str(proj_str)
    delta = float(delta)

    N = lon.size
    if lat.size != N:
        raise RuntimeError("Size of lon and lat must be equal.")

    # Determine the projection
    proj_split = [p.split("=") for p in proj_str.split()]
    strip_plus = lambda s : s[1:] if len(s) > 0 and s[0] == '+' else s
    params = {strip_plus(p[0]) : p[1] for p in proj_split if len(p) > 1}
    if "proj" not in params:
        raise RuntimeError("No projection given.")

    # Create the output buffers and call C code:
    gradient_east = np.zeros((N,2))
    gradient_north = np.zeros((N,2))
    proj_str = proj_str.encode("ascii")
    res = _cdll.gradients_east_north(c_char_p(proj_str), c_ulong(N),
                            lon.ctypes.data_as(POINTER(c_double)),
                            lat.ctypes.data_as(POINTER(c_double)),
                            gradient_east.ctypes.data_as(POINTER(c_double)),
                            gradient_north.ctypes.data_as(POINTER(c_double)),
                            c_double(delta))

    if res != 0:
        raise RuntimeError("grad_east_north backend failed.")

    return gradient_east, gradient_north


def scale_factors(proj_str: str, lon: np.ndarray, lat: np.ndarray,
                  delta: float = 1e-5) -> tuple[np.ndarray,np.ndarray]:
    """

    """
    # Sanity:
    lon = np.ascontiguousarray(lon, dtype=np.double)
    lat = np.ascontiguousarray(lat, dtype=np.double)
    proj_str = str(proj_str)
    delta = float(delta)

    N = lon.size
    if lat.size != N:
        raise RuntimeError("Size of lon and lat must be equal.")

    # Determine the projection
    proj_split = [p.split("=") for p in proj_str.split()]
    strip_plus = lambda s : s[1:] if len(s) > 0 and s[0] == '+' else s
    params = {strip_plus(p[0]) : p[1] for p in proj_split if len(p) > 1}
    if "proj" not in params:
        raise RuntimeError("No projection given.")

    # Create the output buffers and call C code:
    kh = np.zeros((N,2))
    proj_str = proj_str.encode("ascii")
    res = _cdll.scale_factors(c_char_p(proj_str), c_ulong(N),
                            lon.ctypes.data_as(POINTER(c_double)),
                            lat.ctypes.data_as(POINTER(c_double)),
                            kh.ctypes.data_as(POINTER(c_double)),
                            c_double(delta))

    if res != 0:
        raise RuntimeError("grad_east_north backend failed.")

    kh = kh.T

    return kh[0,:], kh[1,:]


