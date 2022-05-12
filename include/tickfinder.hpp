/*
 * Utility for finding ticks along an axis.
 */

#include <../include/types.hpp>
#include <../include/projwrapper.hpp>
#include <../include/griddedinverter.hpp>
#include <vector>

#ifndef PROJPLOT_TICKFINDER_HPP
#define PROJPLOT_TICKFINDER_HPP

namespace projplot {

std::vector<geo_degrees_t>
compute_ticks(const ProjWrapper&, const GriddedInverter&, axis_t, tick_t,
              double xmin, double xmax, double ymin, double ymax,
              double tick_spacing);

}

#endif
