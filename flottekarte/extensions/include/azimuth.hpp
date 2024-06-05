/*
 * A method to disentangle an azimuth field, hopefully yielding a pleasantly
 * smoothed field.
 *
 * Authors: Malte J. Ziebarth (malte.ziebarth@tum.de)
 *
 * Copyright (C) 2024 Technische Universität München
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

#ifndef FLOTTEKARTE_AZIMUTH_HPP
#define FLOTTEKARTE_AZIMUTH_HPP

#include <stddef.h>
#include <cstdint>

namespace flottekarte {

void unwrap_azimuth_field(
    double* angle, uint32_t nx, uint32_t ny, size_t Nmax,
    double cost_beta
);

}

#endif