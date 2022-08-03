/*
 * Compute the boundary of a map based on the (approximate) bijectivity of the
 * map projection.
 *
 * Authors: Malte J. Ziebarth (ziebarth@gfz-potsdam.de)
 *
 * Copyright (C) 2022 Malte J. Ziebarth
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

#include <../include/boundary.hpp>
#include <unordered_map>

using flottekarte::ProjWrapper;
using flottekarte::ProjError;
using flottekarte::path_xy_t;

constexpr double BIJECTIVE_TOLERANCE = 1e-3;
constexpr size_t BOUNDARY_SAMPLES = 100;

#include <iostream>

void flottekarte::bounding_polygon(const ProjWrapper&  proj, double xmin,
                                   double xmax, double ymin, double ymax,
                                   path_xy_t& poly)
{
	/* The method that checks whether a point is part of the map area: */
	auto inside = [&](const xy_t& xy) -> bool {
		try {
			xy_t xy2(proj.project(proj.inverse(xy)));
			const double dx = xy.x - xy2.x;
			const double dy = xy.y - xy2.y;
			return dx*dx + dy*dy <= BIJECTIVE_TOLERANCE * BIJECTIVE_TOLERANCE;
		} catch (const ProjError& err) {
			return false;
		}
	};

	const double Dx = xmax - xmin;
	const double Dy = ymax - ymin;

	/* Generate a point on the bounding box: */
	auto boundary_point = [&](uint_fast8_t boundary, size_t j) -> xy_t {
		switch (boundary){
			case 0:
				/* Bottom */
				return xy_t(xmin + j*Dx/BOUNDARY_SAMPLES, ymin);
			case 1:
				/* Right */
				return xy_t(xmax, ymin + j*Dy/BOUNDARY_SAMPLES);
			case 2:
				/* Top */
				return xy_t(xmax - j*Dx/BOUNDARY_SAMPLES, ymax);
			default:
				/* Left */
				return xy_t(xmin, ymax - j*Dy/BOUNDARY_SAMPLES);
		}
	};

	/* Check whether the whole boundary is inside.
	 * We remember this info since the call to `inside` can be quite expensive.
	 */
	path_xy_t boundary;
	bool all_boundary_in = true;
	for (uint_fast8_t i=0; i<4; ++i){
		for (uint_fast16_t j=0; j<BOUNDARY_SAMPLES; ++j){
			xy_t xy_ij(boundary_point(i,j));
			bool bij = inside(xy_ij);
			all_boundary_in &= bij;
			if (bij)
				boundary.push_back(xy_ij);
		}
	}

	/* Shortcut if all boundary inside: */
	if (all_boundary_in){
		poly.resize(4);
		poly[0] = xy_t(xmin,ymin);
		poly[1] = xy_t(xmax,ymin);
		poly[2] = xy_t(xmax,ymax);
		poly[3] = xy_t(xmin,ymax);
		return;
	}

	/* Check if any boundary is inside: */
	if (!boundary.empty()){
		boundary.swap(poly);
	}

}
