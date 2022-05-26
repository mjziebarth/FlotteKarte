/*
 * Some common types.
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
