#!/bin/python
# This script checks whether the compiled backend `libflottekarte.so`
# can be loaded (or, alternatively, if there is a linking error).
#
# Authors: Malte J. Ziebarth (ziebarth@gfz-potsdam.de)
#
# Copyright (C) 2022 Malte J. Ziebarth
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

import pathlib
from ctypes import CDLL

# Load the shared library:
_parent_directory = pathlib.Path(__file__).parent
_cdll_path = _parent_directory / 'libflottekarte.so'
try:
    CDLL(_cdll_path)
except OSError:
    raise ImportError("Could not load the compiled backend "
                      "'libflottekarte.so'.")
