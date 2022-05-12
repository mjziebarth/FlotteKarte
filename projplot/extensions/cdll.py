# CDLL.

import pathlib
import numpy as np
from typing import Optional
from ctypes import CDLL, c_double, c_size_t, c_uint, POINTER, c_ushort,\
                   c_char_p, c_ulong

# Load the shared library:
_cdll_path = pathlib.Path(__file__).parent / 'libinverse.so'
_cdll = CDLL(_cdll_path)
