/*
 * Path processing facilities.
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

#include <../include/types.hpp>

#ifndef FLOTTEKARTE_PATHS_HPP
#define FLOTTEKARTE_PATHS_HPP

namespace flottekarte {

xy_t boundary_intersection(const xy_t& out, const xy_t& in, double xmin,
                           double xmax, double ymin, double ymax);

struct refined_path_t {
	std::vector<path_xy_t> segments;
	std::vector<xy_t> cuts;
};

refined_path_t
crop_and_refine(const path_geo_t& geo_path, const ProjWrapper& proj,
                double xmin, double xmax, double ymin, double ymax,
                double bisection_offset, double minimum_node_distance);



}

#endif