/*
 * Some common types.
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
#include <projwrapper.hpp>

#ifndef FLOTTEKARTE_TYPES_HPP
#define FLOTTEKARTE_TYPES_HPP

namespace flottekarte {

/* Use defines from projwrappertypes.hpp: */
constexpr double PI = projwrapper::PI;
typedef projwrapper::xy_t xy_t;
typedef projwrapper::geo_t geo_t;
typedef projwrapper::geo_degrees_t geo_degrees_t;
using projwrapper::deg2rad;
using projwrapper::rad2deg;
using projwrapper::modulo;
using projwrapper::ProjWrapper;
using projwrapper::ProjError;


typedef std::vector<geo_degrees_t> path_geo_t;
typedef std::vector<xy_t> path_xy_t;

enum axis_t {
	AX_BOT=0, AX_TOP=1, AX_LEFT=2, AX_RIGHT=3
};

enum tick_t {
	TICK_LON=0, TICK_LAT=1, TICK_NONE=2
};

}

#endif
