/*
 * Utility for finding ticks along an axis.
 */

#include <../include/types.hpp>
#include <../include/projwrapper.hpp>
#include <../include/griddedinverter.hpp>
#include <vector>

#ifndef FLOTTEKARTE_TICKFINDER_HPP
#define FLOTTEKARTE_TICKFINDER_HPP

namespace flottekarte {

std::vector<geo_degrees_t>
compute_ticks(const ProjWrapper&, const GriddedInverter&, axis_t, tick_t,
              double xmin, double xmax, double ymin, double ymax,
              double tick_spacing);

}

#endif
