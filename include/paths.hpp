/*
 * Path processing facilities.
 */

#include <../include/types.hpp>
#include <../include/projwrapper.hpp>

#ifndef PROJPLOT_PATHS_HPP
#define PROJPLOT_PATHS_HPP

namespace projplot {

std::vector<path_xy_t>
crop_and_refine(const path_geo_t& geo_path, const ProjWrapper& proj,
                double xmin, double xmax, double ymin, double ymax,
                double bisection_offset, double minimum_node_distance,
                double cut_at_angle_degrees, bool cyclic);

}

#endif