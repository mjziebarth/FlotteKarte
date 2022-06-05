# Load GeoJSON data.

import json
import numpy as np
from pyproj import Proj
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

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
    def __init__(self, fname, proj):
        # Initialize the projection:
        if not isinstance(proj,Proj):
            proj = Proj(proj)

        # Load all content:
        with open(fname,'r') as f:
            geojson = json.load(f)

        # Some sanity assertion:
        assert_geojson_sanity(geojson)

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
                        poly.append(np.stack((x,y),axis=1))
                    multipoly.append(poly)
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


    def get_polygon_patches(self):
        """
        Get matplotlib patches of the polygons and multipolygons.
        """
        poly_patches = []
        for poly in self.polygons:
            poly_patches.append(Polygon(poly[0]))
        for mpoly in self.multipolygons:
            for poly in mpoly:
                poly_patches.append(Polygon(poly[0]))

        return PatchCollection(poly_patches)


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
