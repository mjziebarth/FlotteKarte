/*
 * Generate a set of paths for a coordinate grid.
 */

#include <vector>
//#include <../include/types.hpp>
#include <../include/projwrapper.hpp>

#ifndef PROJPLOT_GRID_HPP
#define PROJPLOT_GRID_HPP

namespace projplot {


std::vector<std::vector<xy_t>>
generate_grid_lines(const ProjWrapper& proj, double xmin, double xmax,
                    double ymin, double ymax, int tick_spacing_degree,
                    double bisection_offset, double minimum_node_distance,
                    double max_lat, double cut_at_angle_degrees);


}

#endif