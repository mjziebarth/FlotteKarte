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

#include <vector>
#include <type_traits>
#include <projwrapper.hpp>

#ifndef FLOTTEKARTE_TYPES_HPP
#define FLOTTEKARTE_TYPES_HPP

namespace flottekarte {

/* Use defines from projwrappertypes.hpp: */
constexpr double PI = projwrapper::PI;
typedef projwrapper::geo_t geo_t;
typedef projwrapper::geo_degrees_t geo_degrees_t;
using projwrapper::deg2rad;
using projwrapper::rad2deg;
using projwrapper::modulo;
using projwrapper::ProjWrapper;
using projwrapper::ProjError;

class xy_t : public projwrapper::xy_t
{
public:
	xy_t(double x, double y) : projwrapper::xy_t(x, y)
	{}

	xy_t() = default;
	xy_t(xy_t&&) = default;
	xy_t(const xy_t&) = default;

	xy_t(projwrapper::xy_t&& xy);
	xy_t(const projwrapper::xy_t& xy);

	xy_t& operator=(const xy_t& other) = default;

	template<
		typename number,
		typename std::enable_if<std::is_arithmetic_v<number>,int>::type = 0
	>
	xy_t operator*(number n) const
	{
		return res(n*x, n*y);
	}


	template<
		typename number,
		typename std::enable_if<std::is_arithmetic_v<number>,int>::type = 0
	>
	xy_t& operator*=(number n)
	{
		x *= n;
		y *= n;
		return *this;
	}

	xy_t operator+(xy_t&& other) const;
	xy_t operator+(const xy_t& other) const;

	xy_t operator-(const xy_t& other) const;

	xy_t& operator+=(const xy_t&);

	double dot(const xy_t& other) const;
};

typedef std::vector<geo_degrees_t> path_geo_t;
typedef std::vector<xy_t> path_xy_t;

enum axis_t {
	AX_BOT=0, AX_TOP=1, AX_LEFT=2, AX_RIGHT=3
};

enum tick_t {
	TICK_LON=0, TICK_LAT=1, TICK_NONE=2
};




/*
 * Further arithmetic operators:
 */
template<
	typename number,
	typename std::enable_if<std::is_arithmetic_v<number>,int>::type = 0
>
xy_t operator*(number n, xy_t&& xy)
{
	xy.x *= n;
	xy.y *= n;
	return xy;
}

template<
	typename number,
	typename std::enable_if<std::is_arithmetic_v<number>,int>::type = 0
>
xy_t operator*(number n, const xy_t& xy)
{
	return xy_t(n*xy.x, n*xy.y);
}

}

#endif
