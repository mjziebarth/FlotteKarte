/*
 * Inversion routine.
 */

#include <../include/projwrapper.hpp>
#include <../include/types.hpp>

#ifndef FLOTTEKARTE_INVERT_HPP
#define FLOTTEKARTE_INVERT_HPP

namespace flottekarte {

geo_t gradient_descent_inverse_project(const ProjWrapper& projection,
                                       const xy_t& xy, double lambda0,
                                       double phi0);

}

#endif
