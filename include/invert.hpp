/*
 * Inversion routine.
 */

#include <../include/projwrapper.hpp>
#include <../include/types.hpp>

#ifndef PROJPLOT_INVERT_HPP
#define PROJPLOT_INVERT_HPP

namespace projplot {

geo_t gradient_descent_inverse_project(const ProjWrapper& projection,
                                       const xy_t& xy, double lambda0,
                                       double phi0);

}

#endif
