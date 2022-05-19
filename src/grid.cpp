/*
 * Generate a set of paths for a coordinate grid.
 */

#include <../include/grid.hpp>
#include <../include/paths.hpp>
#include <cmath>
#include <iostream>

using projplot::xy_t;
using projplot::geo_degrees_t;
using projplot::generate_grid_lines;
using projplot::ProjWrapper;
using projplot::crop_and_refine;
using projplot::path_xy_t;
using projplot::path_geo_t;


static long lon_ticks(int tick_spacing_degree)
{
	if (tick_spacing_degree > 0){
		return std::max<long>(std::floor(180.0 / tick_spacing_degree), 1);
	} else {
		return 180 * static_cast<long>(-tick_spacing_degree);
	}
}

static long lat_ticks(int tick_spacing_degree)
{
	if (tick_spacing_degree > 0){
		return std::max<size_t>(std::floor(90.0 / tick_spacing_degree), 1);
	} else {
		return 90 * static_cast<size_t>(-tick_spacing_degree);
	}
}



std::vector<path_xy_t>
projplot::generate_grid_lines(const ProjWrapper& proj, double xmin, double xmax,
                    double ymin, double ymax, int tick_spacing_degree,
                    double bisection_offset, double minimum_node_distance,
                    double max_lat, double cut_at_angle_degrees)
{
	if (tick_spacing_degree == 0){
		throw std::runtime_error("Cannot have zero ticks spacing.");
	}
	if (max_lat <= 0){
		throw std::runtime_error("Cannot have negative or zero max_lat.");
	}
	max_lat = std::min(max_lat, 90.0);
	const double diag = std::sqrt(  (xmax-xmin)*(xmax-xmin)
	                              + (ymax-ymin)*(ymax-ymin));
	/* Want a minimum number of 10 ticks per diagonal area of the map
	 * extent: */
	const size_t min_ticks_per_degree
	    = std::ceil(PI/180.0*proj.a() / diag * 10);
	if (min_ticks_per_degree > 1e5){
		/* The map is too small so more than 1e5 ticks per degree are
		 * required to have ~10 ticks within the map area (per line).
		 * Here, we would have to use a more efficient algorithm that
		 * takes into account which geographic coordinates can occur
		 * within the map extents. */
		throw std::runtime_error("Too many ticks requested. Abort!");
	}


	/* Determine the number of ticks in longitude direction from 0 to 180°: */
	const long nlon0 = lon_ticks(tick_spacing_degree);
	const double dlon0 = (tick_spacing_degree > 0) ? tick_spacing_degree
	                          : 1.0 / tick_spacing_degree;

	/* Determine the number of ticks in latitude direction from 0 to 90°: */
	const long nlat1 = lat_ticks(tick_spacing_degree);
	const double dlat1 = (tick_spacing_degree > 0) ? tick_spacing_degree
	                          : 1.0 / tick_spacing_degree;


	std::vector<path_xy_t> paths;
	/* Construct the meridian lines: */
	const long nlat0 = std::ceil(max_lat*min_ticks_per_degree);
	const double dlat0 = max_lat / nlat0;
	path_geo_t path;
	for (long i = -nlon0+1; i<nlon0; ++i){
		/* Construct the vector of geographic coordinates
		 * from -max_lat to max_lat:
		 */
		path.clear();
		const double loni = i*dlon0;
		for (long j=-nlat0+1; j<nlat0; ++j){
			path.emplace_back(loni, j*dlat0);
		}

		/* Crop and refine the path: */
		std::vector<path_xy_t>
		   refined(crop_and_refine(path, proj, xmin, xmax, ymin, ymax,
		                           bisection_offset, minimum_node_distance,
		                           cut_at_angle_degrees, false));

		/* Add the new segments: */
		paths.insert(paths.cend(), refined.cbegin(), refined.cend());
	}

	/* Construct the lines of constant latitude: */
	const long nlon1 = 180 * min_ticks_per_degree;
	const double dlon1 = 180.0 / nlon1;
	path.clear();
	for (long i = -nlat1+1; i<nlat1; ++i){
		/* Construct the vector of geographic coordinates
		 * from -max_lat to max_lat:
		 */
		path.clear();
		const double lati = i*dlat1;
		for (long j=-nlon1; j<nlon1; ++j){
			path.emplace_back(j*dlon1, lati);
		}

		/* Crop and refine the path: */
		std::vector<path_xy_t>
		   refined(crop_and_refine(path, proj, xmin, xmax, ymin, ymax,
		                           bisection_offset, minimum_node_distance,
		                           cut_at_angle_degrees, true));

		/* Add the new segments: */
		paths.insert(paths.cend(), refined.cbegin(), refined.cend());
	}

	/* End: */
	return paths;
}
