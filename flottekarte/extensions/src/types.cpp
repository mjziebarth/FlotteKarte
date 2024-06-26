/*
 * Some common types.
 *
 * Authors: Malte J. Ziebarth (ziebarth@gfz-potsdam.de)
 *
 * Copyright (C) 2022 Deutsches GeoForschungsZentrum Potsdam,
 *               2024 Technische Universität München
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
#include <cmath>

namespace flottekarte {

xy_t::xy_t(projwrapper::xy_t&& xy) : projwrapper::xy_t(std::move(xy))
{}

xy_t::xy_t(const projwrapper::xy_t& xy) : projwrapper::xy_t(xy)
{}

xy_t& xy_t::operator+=(const xy_t& other)
{
    x += other.x;
    y += other.y;
    return *this;
}

xy_t xy_t::operator+(xy_t&& other) const
{
    other.x += x;
    other.y += y;
    return other;
}

xy_t xy_t::operator+(const xy_t& other) const
{
    return xy_t(x + other.x, y + other.y);
}

xy_t xy_t::operator-(const xy_t& other) const
{
    return xy_t(x - other.x, y - other.y);
}

double xy_t::dot(const xy_t& other) const
{
    return x * other.x + y * other.y;
}

double xy_t::norm() const
{
    return std::sqrt(x*x + y*y);
}

} // end namespace
