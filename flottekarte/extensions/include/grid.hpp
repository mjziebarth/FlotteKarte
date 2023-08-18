/*
 * Generate a set of paths for a coordinate grid.
 *
 * Authors: Malte J. Ziebarth (ziebarth@gfz-potsdam.de)
 *
 * Copyright (C) 2022 Deutsches GeoForschungsZentrum Potsdam
 *
 * Licensed under the EUPL, Version 1.2 or – as soon they will be approved by
 * the European Commission - subsequent versions of the EUPL (the "Licence");
 * You may not use this work except in compliance with the Licence.
 * You may obtain a copy of the Licence at:
 *
 * https://joinup.ec.europa.eu/collection/eupl/eupl-text-eupl-12
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the Licence is distributed on an "AS IS" basis,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the Licence for the specific language governing permissions and
 * limitations under the Licence.
 */

#include <vector>
#include <../include/types.hpp>

#ifndef FLOTTEKARTE_GRID_HPP
#define FLOTTEKARTE_GRID_HPP

namespace flottekarte {

struct cut_t {
	xy_t point;
	tick_t tick_type;
	double coordinate;
};

struct geo_grid_t {
	std::vector<path_xy_t> paths;
	std::vector<cut_t> cuts;
};


geo_grid_t
generate_grid_lines(const ProjWrapper& proj, double xmin, double xmax,
                    double ymin, double ymax, int tick_spacing_degree,
                    double bisection_offset, double minimum_node_distance,
                    double max_lat);


}

#endif