# Testing the GeoJSON data loading facility.
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
from flottekarte import GeoJSON
from pathlib import Path

def test_load_iceland():
    """
    Tries to load the Natural Earth iceland polygon
    (data/natural-earth/iceland.geojson).
    """
    path = Path(__file__).parent / 'data/natural-earth/iceland.geojson'
    geojson = GeoJSON(path.resolve(), '+proj=stere +lon_0=-18 lat_0=64.5')

    # Some previously computed attributes:
    assert len(geojson.multipolygons) == 1
    assert len(geojson.multipolygons[0]) == 5
    assert np.allclose(geojson.get_extent(),
                       [(-302210.1135546372, 211360.90840587052),
                        (-122781.37705651142, 230189.97898052254)])

    # The rectangle defined by xlim=(50e3, 100e3), ylim=(100e3, 200e3)
    # is intercepted only by the main Iceland shoreline and does not
    # contain any of the other 4 islands.
    # Hence, we can test the filtering capability here:
    geojson2 = GeoJSON(path.resolve(), '+proj=stere +lon_0=-18 lat_0=64.5',
                       xlim=(50e3, 100e3),
                       ylim=(100e3, 200e3))
    assert len(geojson2.multipolygons) == 1
    assert len(geojson2.multipolygons[0]) == 1
