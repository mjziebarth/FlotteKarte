/*
 * Coordinate inversion based on gradient descent of mismatch to PROJ forward
 * projection.
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

#ifndef FLOTTEKARTE_INVERT_HPP
#define FLOTTEKARTE_INVERT_HPP

namespace flottekarte {

geo_t gradient_descent_inverse_project(const ProjWrapper& projection,
                                       const xy_t& xy, double lambda0,
                                       double phi0);

}

#endif
