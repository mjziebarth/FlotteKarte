# Load GeoJSON data.
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

import json
import numpy as np
from pyproj import Proj
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from warnings import warn

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
    """
    def __init__(self, fname, proj, xlim=None, ylim=None):
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
                    poly.append(np.stack((x,y),axis=1))
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
                        poly.append(np.stack((x,y),axis=1))
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
