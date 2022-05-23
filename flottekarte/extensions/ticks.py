# Axes ticks computation.

import numpy as np
from .cdll import _cdll
from ctypes import c_double, POINTER, c_char_p, c_uint, c_ubyte
from typing import Optional

def compute_axes_ticks(proj_str: str, xmin: float, xmax: float, ymin: float,
                       ymax: float, tick_spacing: float,
                       bot: Optional[str] = 'lon', top: Optional[str] = 'lon',
                       left: Optional[str] = 'lat',
                       right: Optional[str] = 'lat') \
   -> tuple[np.ndarray,np.ndarray,np.ndarray,np.ndarray]:
    """

    """
    # Sanity:
    xmin = float(xmin)
    xmax = float(xmax)
    ymin = float(ymin)
    ymax = float(ymax)
    tick_spacing = float(tick_spacing)
    proj_str = str(proj_str)
    t2i = {'lon' : 0, 'lat': 1, None: 2}
    if any(t not in t2i for t in (bot,top,left,right)):
        raise RuntimeError("Ticks `bot`, `top`, `left` and `right` must "
                           "be one of 'lon', 'lat', or None.")
    bot = t2i[bot]
    top = t2i[top]
    left = t2i[left]
    right = t2i[right]

    # Determine the projection
    proj_split = [p.split("=") for p in proj_str.split()]
    strip_plus = lambda s : s[1:] if len(s) > 0 and s[0] == '+' else s
    params = {strip_plus(p[0]) : p[1] for p in proj_split if len(p) > 1}
    if "proj" not in params:
        raise RuntimeError("No projection given.")

    # Create the output buffers and call C code:
    Nmax = 100
    ticks_bot = np.zeros((Nmax,2))
    ticks_top = np.zeros((Nmax,2))
    ticks_left = np.zeros((Nmax,2))
    ticks_right = np.zeros((Nmax,2))
    proj_str = proj_str.encode("ascii")
    ticks_lengths = np.zeros(4, dtype=np.uintc)

    res = _cdll.compute_axes_ticks(c_char_p(proj_str), c_double(xmin),
                            c_double(xmax), c_double(ymin), c_double(ymax),
                            c_double(tick_spacing), c_uint(Nmax), c_ubyte(bot),
                            c_ubyte(top), c_ubyte(left), c_ubyte(right),
                            ticks_bot.ctypes.data_as(POINTER(c_double)),
                            ticks_top.ctypes.data_as(POINTER(c_double)),
                            ticks_left.ctypes.data_as(POINTER(c_double)),
                            ticks_right.ctypes.data_as(POINTER(c_double)),
                            ticks_lengths.ctypes.data_as(POINTER(c_uint)))

    if res != 0:
        raise RuntimeError("grad_east_north backend failed.")

    return (ticks_bot[:ticks_lengths[0]], ticks_top[:ticks_lengths[1]],
            ticks_left[:ticks_lengths[2]], ticks_right[:ticks_lengths[3]])
