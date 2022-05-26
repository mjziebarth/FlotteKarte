/*
 * Generate a set of paths for a coordinate grid.
 */

#include <vector>
#include <../include/types.hpp>

#ifndef FLOTTEKARTE_GRID_HPP
#define FLOTTEKARTE_GRID_HPP

namespace flottekarte {


std::vector<std::vector<xy_t>>
generate_grid_lines(const ProjWrapper& proj, double xmin, double xmax,
                    double ymin, double ymax, int tick_spacing_degree,
                    double bisection_offset, double minimum_node_distance,
                    double max_lat, double cut_at_angle_degrees);


}

#endif