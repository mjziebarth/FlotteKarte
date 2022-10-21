# Computing automatic bounds for map projections.
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

from typing import Tuple, Optional
from pyproj import Proj
import numpy as np

_GLOBAL_PROJECTIONS_1 = set(["adams_ws1", "adams_ws2", "aitoff", "bacon",
                             "boggs", "cea", "comill", "crast", "denoy", "eck1",
                             "eck2", "eck3", "eck4", "eck5", "eck6", "eqearth",
                             "fahey", "fouc", "fouc_s", "gall", "gins8",
                             "gn_sinu", "goode", "hammer", "hatano", "igh",
                             "igh_o", "kav5", "kav7", "lagrng", "loxim",
                             "mbt_s", "mbt_fps", "mbtfpp", "mbtfpq", "mbtfps",
                             "mill", "moll", "natearth", "natearth2", "nell",
                             "nell_h", "nicol", "ortel", "patterson", "putp1",
                             "putp2", "putp3", "putp3p", "putp4p", "putp5",
                             "putp5p", "putp6", "putp6p", "qua_aut", "robin",
                             "sinu", "times", "urm5", "vandg", "vandg2",
                             "vandg3", "wag1", "wag2", "wag3", "wag4", "wag5",
                             "wag6", "weren", "wink1", "wink2", "winktri"])


def automatic_map_extents(proj_str: str, proj: Optional[Proj] = None)\
    -> Tuple[float,float,float,float]:
    """
    Try to compute map extents automatically.
    """
    # Compute the PROJ instance:
    if proj is None:
        proj = Proj(proj_str)

    # Extracting parameters from a proj string:
    def extract_proj_parameter(proj_str: str, parameter: str) -> str:
        """
        Extracts a proj parameter.
        """
        p = parameter + "="
        return [s for s in proj_str.split() if p in s][0].split(p)[1]


    # Extract the map projection definition:
    projection = extract_proj_parameter(proj_str, "proj")

    # Now handle the known cases:
    if projection in _GLOBAL_PROJECTIONS_1:
        if "lon_0" in proj_str:
            lon_0 = float(extract_proj_parameter(proj_str, "lon_0"))
        else:
            lon_0 = 0.0
        xmin = proj((lon_0-179.999) % 360.0, 0)[0]
        xmax = proj((lon_0+179.999) % 360.0, 0)[0]
        ymin = proj(lon_0, -90.0)[1]
        ymax = proj(lon_0, 90.0)[1]

    elif projection in ["larr","wag7"]:
        lon_0 = float(extract_proj_parameter(proj_str, "lon_0"))
        # Evaluate at +/- 180° longitude relative to lon_0.
        lat = np.concatenate([np.linspace(-90,-80,50), [0],
                              np.linspace(80,90,50)])
        x,y = proj(np.ones(101) * (lon_0-179.999) % 360.0, lat)
        xmin = x.min()
        ymin = y.min()
        ymax = y.max()
        x,y = proj(np.ones(101) * (lon_0+179.999) % 360.0, lat)
        xmax = x.max()
        ymin = min(ymin,y.min())
        ymax = max(ymax,y.max())

    elif projection == "collg":
        lon_0 = float(extract_proj_parameter(proj_str, "lon_0"))
        xmin,ymin = proj((lon_0-179.999) % 360.0, -90)
        xmax = proj((lon_0+179.999) % 360.0, -90)[0]
        ymax = proj(lon_0, 90.0)[1]

    else:
        raise NotImplementedError("Cannot determine map extent for projection\""
                                  + str(proj_str) + "\" automatically.")

    return xmin, xmax, ymin, ymax
