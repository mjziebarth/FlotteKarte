# CDLL.

import pathlib
import numpy as np
from ctypes import CDLL

# Load the shared library:
_cdll_path = pathlib.Path(__file__).parent / 'libinverse.so'
_cdll = CDLL(_cdll_path)
