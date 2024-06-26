/*
 * Two-dimensional interpolation.
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
 * */

#include <../include/interpolate.hpp>

namespace flottekarte {

grid2d_t<true>::grid2d_t(
    double xmin,
    double xmax,
    size_t nx,
    double ymin,
    double ymax,
    size_t ny,
    const double* z
)
   : xmin(xmin), xmax(xmax), ymin(ymin), ymax(ymax), nx(nx), ny(ny),
     data(z, z+nx*ny), z(data.data())
{}

grid2d_t<false>::grid2d_t(
    double xmin,
    double xmax,
    size_t nx,
    double ymin,
    double ymax,
    size_t ny,
    const double* z
)
   : xmin(xmin), xmax(xmax), ymin(ymin), ymax(ymax), nx(nx), ny(ny),
     z(z)
{}







}