# CDLL-based interfacing with the C++ code.
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

import os
import pathlib
import subprocess
import numpy as np
from ctypes import CDLL
from shutil import copyfile
from warnings import warn

# Load the shared library:
_parent_directory = pathlib.Path(__file__).parent
_cdll_path = _parent_directory / 'libflottekarte.so'
try:
    _cdll = CDLL(_cdll_path)
except OSError:
    _cdll = None
    raise ImportError("Could not load the compiled backend "
                      "'libflottekarte.so'. Most likely this means that your "
                      "system library has been updated and you need to "
                      "recompile the FlotteKarte backend. You can do so from "
                      "within Python by importing and calling the function "
                      "`flottekarte.recompile_flottekarte`. Afterwards, "
                      "restarting Python will be necessary.")


# Recompilation facility:
def recompile_flottekarte():
    """
    Recompile the C++ backend and link against updated PROJ system
    library.
    """
    global _cdll
    if _cdll is not None:
        warn("Recompiling with a loaded backend can lead to a Python crash. "
             "Please restart Python after the compilation was successful.")
    else:
        _cdll = None

    # Current directory:
    current_dir = pathlib.Path.cwd().absolute()

    # Path to this file. In the parent directory, we should find copies
    # of the `include` and `src` folders.
    include = _parent_directory / "include"
    src = _parent_directory / "src"
    subproj = _parent_directory / "subprojects"
    if not include.is_dir() or not src.is_dir() or not subproj.is_dir():
        raise RuntimeError("One of the source directories has not been "
                           "copied to the extension path. Cannot compile.")
    if not (_parent_directory / "meson.build").is_file():
        raise RuntimeError("Did not find `meson.build` file in extension "
                           "path.")

    # Change to the extension path:
    os.chdir(_parent_directory)

    # Perform Meson build:
    subprocess.run(["meson","setup","builddir"], check=True)
    subprocess.run(["meson","compile","-C","builddir"], check=True)

    # Copy the compiled library:
    copyfile((_parent_directory / "builddir" / "libflottekarte.so").absolute(),
             (_parent_directory / "libflottekarte.so").absolute())

    # Ensure that we do not try to access any loaded CDLL:
    _cdll = None

    # Change back to working directory:
    os.chdir(current_dir)
