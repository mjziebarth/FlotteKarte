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
from shutil import copyfile, rmtree
from logging import info

# Paths to the extension:
_parent_directory = pathlib.Path(__file__).parent
_cdll_path = _parent_directory / 'libflottekarte.so'


#
# Info for system packaging this Python package:
# This file would have to be modified since recompilation should not write
# into the system-wide Python packages directory.
# The meson build would have to be performed by the packager, and
# libflottekarte.so would have to be shipped with the package.
# The following code can all be removed within [BEGIN] to [END]
#
# [BEGIN]

# Whether to keep the Meson build directory or not.
_keep_meson_builddir_subprojects = False


# A function to check whether the flottekarte shared object can
# be loaded. Does so in a separate process so that, on failure,
# the object can be replaced:
def can_load_flottekarte():
    """
    This function checks whether the compiled backend can be loaded.
    """
    try:
        res = subprocess.run(["python",_parent_directory / "test-cdll.py"],
                             check=True, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError:
        return False
    return True


# Recompilation facility:
def recompile_flottekarte(verbose: bool = False):
    """
    Recompile the C++ backend and link against updated PROJ system
    library.
    """
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
    kwargs = {}
    if not verbose:
        kwargs["stdout"] = subprocess.DEVNULL
        kwargs["stderr"] = subprocess.DEVNULL
    subprocess.run(["meson","setup","builddir"], check=True, **kwargs)
    subprocess.run(["meson","compile","-C","builddir"], check=True, **kwargs)

    # Copy the compiled library:
    copyfile((_parent_directory / "builddir" / "libflottekarte.so").absolute(),
             (_parent_directory / "libflottekarte.so").absolute())

    if not _keep_meson_builddir_subprojects:
        rmtree((_parent_directory / "builddir").absolute())
        rmtree((_parent_directory / "subprojects" / "libprojwrap").absolute())

    # Change back to working directory:
    os.chdir(current_dir)


# See if the import works:
if not can_load_flottekarte():
    # Recompile:
    print("Need to recompile FlotteKarte backend libflottekarte.so...")
    recompile_flottekarte()

# [END]

# Load the shared library:
try:
    _cdll = CDLL(_cdll_path)
except OSError:
    raise ImportError("Could not load the compiled backend "
                      "'libflottekarte.so'. Trying to recompile failed.")
