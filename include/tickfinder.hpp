/*
 * Utility for finding ticks along an axis.
 *
 * Authors: Malte J. Ziebarth (ziebarth@gfz-potsdam.de)
 *
 * Copyright (C) 2022 Deutsches GeoForschungsZentrum Potsdam,
 *                    Malte J. Ziebarth
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
#include <../include/augmentedproj.hpp>
#include <vector>

#ifndef FLOTTEKARTE_TICKFINDER_HPP
#define FLOTTEKARTE_TICKFINDER_HPP

namespace flottekarte {

struct segment_tick_t {
	size_t segment;
	geo_degrees_t tick;

	segment_tick_t() = default;

	segment_tick_t(size_t segment, double lon, double lat);
};

std::vector<segment_tick_t>
compute_ticks(const AugmentedProj&, tick_t,
              const path_xy_t&, double tick_spacing);

}

#endif
