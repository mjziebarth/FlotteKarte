# Load GeoJSON data.
#
# Authors: Malte J. Ziebarth (ziebarth@gfz-potsdam.de)
#
# Copyright (C) 2022 Deutsches GeoForschungsZentrum Potsdam,
#               2022 Malte J. Ziebarth
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

import json
import numpy as np
from typing import Union, Optional, Tuple
from pyproj import Proj
from matplotlib.path import Path
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from warnings import warn

# Using matplotlib backend to clip the polygon paths to the map extent:
try:
    from matplotlib import _path
    _has_matplotlib_path = True
except ImportError:
    _has_matplotlib_path = False
    warn("Could not load the matplotlib '_path' module. Polygon clipping "
         "to map extents will not be possible.")


def assert_geojson_sanity(geojson_dict):
    """
    This function tests some basic sanity conditions of the
    GeoJSON.

    Note: Does not necessarily follow all possibilities of
    the standard (yet).
    """
    assert geojson_dict["type"] == "FeatureCollection"
    assert "features" in geojson_dict


class GeoJSON:
    """
    Loads data from a GeoJSON.

    Parameters
    ----------
    fname : str
       Path to the GeoJSON file. Alternatively
       anything that can be accepted by :func:`open`.
    proj : str or pyproj.Proj
       Projection to evaluate the geometries in.
    xlim : tuple[float,float], optional
       Selection *x* extent in projected coordinates.
       Given as :math:`(x_\\mathrm{min},x_\\mathrm{max})`.
    ylim : tuple[float,float], optional
       Selection *y* extent in projected coordinates.
       Given as :math:`(y_\\mathrm{min},y_\\mathrm{max})`.


    Notes
    -----
    If no extents are given, all geometries are loaded.
    So far, only polygon and multipolygon geometries are
    loaded.
    """
    def __init__(self, fname: str, proj: Union[str,Proj],
                 xlim: Optional[Tuple[float,float]] = None,
                 ylim: Optional[Tuple[float,float]] = None):
        # Initialize the projection:
        if not isinstance(proj,Proj):
            proj = Proj(proj)

        # Load all content:
        with open(fname,'r') as f:
            geojson = json.load(f)

        # Some sanity assertion:
        assert_geojson_sanity(geojson)

        # Simple extent selection:
        if xlim is not None:
            if ylim is not None:
                def select_points(x,y):
                    return np.any((x >= xlim[0]) & (x <= xlim[1]) &
                                  (y >= ylim[0]) & (y <= ylim[1]))
            else:
                def select_points(x,y):
                    return np.any((x >= xlim[0]) & (x <= xlim[1]))
            # TODO FIXME.
            warn("In GeoJSON, polygon selection is currently point selection. "
                 "This can cause polygons not to be selected if they "
                 "completely contain the x and y extents.")
            select_polygon = select_points
        elif ylim is not None:
            def select_points(x,y):
                    return np.any((y >= ylim[0]) & (y <= ylim[1]))
            # TODO FIXME.
            warn("In GeoJSON, polygon selection is currently point selection. "
                 "This can cause polygons not to be selected if they "
                 "completely contain the x and y extents.")
            select_polygon = select_points
        else:
            def select_points(x,y):
                return True
            select_polygon = select_points

        # Cropping the paths:
        if _has_matplotlib_path and xlim is not None and ylim is not None:
            def crop_poly(poly):
                # Create an (N+1,2)-shaped array of coordinates,
                # with the last entry being unused (corresponds to
                # matplotlib Path.CLOSEPOLY code)
                path = np.empty((poly.shape[0]+1, 2))
                path[:-1,:] = poly

                # Code of a closed matplotlib Path:
                codes = np.full(path.shape[0], Path.LINETO,
                                dtype=Path.code_type)
                codes[0] = Path.MOVETO
                codes[-1] = Path.CLOSEPOLY
                path = Path(path, codes=codes)

                # Use the Matplotlib C++ backend code to crop the path
                # to the map extents:
                return _path.clip_path_to_rect(path, ((xlim[0],ylim[0]),
                                                      (xlim[1],ylim[1])),
                                               True)
        else:
            # No-op:
            def crop_poly(poly):
                return [poly]

        # Iterate through the features:
        points = []
        multipoints = []
        linestrings = []
        multilinestrings = []
        polygons = []
        multipolygons = []
        for feat in geojson["features"]:
            assert feat["type"] == "Feature"
            geom = feat["geometry"]
            geom_type = geom["type"]
            if geom_type ==  "Point":
                pass
            elif geom_type == "MultiPoint":
                pass
            elif geom_type == "LineString":
                pass
            elif geom_type == "MultiLineString":
                pass
            elif geom_type == "Polygon":
                # Read and project a polygon.
                rings = geom["coordinates"]
                poly = []
                for lola in rings:
                    x,y = proj(*np.array(lola).T)
                    if not select_polygon(x,y):
                        continue
                    poly.extend(crop_poly(np.stack((x,y), axis=1)))
                polygons.append(poly)

            elif geom_type == "MultiPolygon":
                # Read and project the polygons.
                feat_polys = geom["coordinates"]
                multipoly = []
                for rings in feat_polys:
                    poly = []
                    for lola in rings:
                        x,y = proj(*np.array(lola).T)
                        if not select_polygon(x,y):
                            continue
                        poly.extend(crop_poly(np.stack((x,y), axis=1)))
                    if len(poly) > 0:
                        multipoly.append(poly)
                if len(multipoly) > 0:
                    multipolygons.append(multipoly)
            else:
                raise RuntimeError("Unknown geometry type detected in "
                                   "GeoJSON feature.")

        self.points = NotImplemented
        self.multipoints = NotImplemented
        self.linestrings = NotImplemented
        self.multilinestrings = NotImplemented
        self.polygons = polygons
        self.multipolygons = multipolygons


    def get_polygon_patches(self, **kwargs):
        """
        Get matplotlib patches of the polygons and multipolygons.

        Parameters
        ----------
        kwargs : optional
           Parameters to pass to the patch collection.

        Returns
        -------
        :class:`matplotlib.collections.PatchCollection`
           Patches for all polygons and multipolygons.
        """
        poly_patches = []
        for poly in self.polygons:
            poly_patches.append(Polygon(poly[0]))
        for mpoly in self.multipolygons:
            for poly in mpoly:
                poly_patches.append(Polygon(poly[0]))

        return PatchCollection(poly_patches, **kwargs)


    def get_extent(self):
        """
        Get x and y limits.

        Returns
        -------
        xlim : tuple[float,float]
           *x* extents :math:`(x_\\mathrm{min},x_\\mathrm{max})`
        ylim : tuple[float,float]
           *y* extents :math:`(y_\\mathrm{min},y_\\mathrm{max})`
        """
        xmin = np.inf
        xmax = -np.inf
        ymin = np.inf
        ymax = -np.inf
        for mpoly in self.multipolygons:
            for poly in mpoly:
                for ring in poly:
                    xmin = min(xmin, ring[:,0].min())
                    xmax = max(xmax, ring[:,0].max())
                    ymin = min(ymin, ring[:,1].min())
                    ymax = max(ymax, ring[:,1].max())

        return (xmin,xmax), (ymin,ymax)
