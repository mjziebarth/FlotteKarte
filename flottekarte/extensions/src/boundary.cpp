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
#include <../include/paths.hpp>
#include <cmath>
#include <stack>

using flottekarte::xy_t;
using flottekarte::AugmentedProj;
using flottekarte::ProjError;
using flottekarte::path_xy_t;
using flottekarte::boundary_intersection;

constexpr double BIJECTIVE_TOLERANCE = 1e-8;
constexpr size_t BOUNDARY_SAMPLES = 101;

#include <iostream>

static inline bool inside(const xy_t& xy, const AugmentedProj& proj)
{
	/* The method that checks whether a point is part of the map area: */
	constexpr double tol = BIJECTIVE_TOLERANCE * BIJECTIVE_TOLERANCE;
	try {
		xy_t xy2(proj.project(proj.inverse(xy)));
		const double dx = xy.x - xy2.x;
		const double dy = xy.y - xy2.y;
		const double norm2 = std::max(xy.x * xy.x + xy.y * xy.y, 1.0);
		return dx*dx + dy*dy <= tol * norm2;
	} catch (const ProjError& err) {
		return false;
	}
}

static xy_t find_boundary(const xy_t& start, double outward_x, double outward_y,
                          const AugmentedProj& proj, double atol)
{
	/* Find the correct search direction from the start: */
	const bool outward = inside(start, proj);
	double search_x, search_y;
	if (outward){
		search_x = outward_x;
		search_y = outward_y;
	} else {
		search_x = -outward_x;
		search_y = -outward_y;
	}

	/* Norm the search vector: */
	const double scale = std::sqrt(search_x*search_x + search_y*search_y);
	if (scale == 0){
		throw std::runtime_error("Encountered zero length segment.");
	}
	search_x /= scale;
	search_y /= scale;

	/* Now define a coordinate `z` from `start` along the direction `search`: */
	auto has_transition = [&](double z) -> bool {
		return inside(xy_t(start.x + search_x * z,
		                   start.y + search_y * z),
		              proj) != outward;
	};

	/* Find the interval that contains the root transition: */
	double zr = atol;
	while (!has_transition(zr)){
		zr *= 2.0;
		if (std::isinf(zr)){
			throw std::runtime_error("Did not find the map boundary.");
		}
	}

	/* Find the root to absolute precision `atol` using bisection: */
	double zl = 0;
	while (zr - zl > atol){
		const double z = 0.5 * (zr + zl);
		if (has_transition(z)){
			zr = z;
		} else {
			zl = z;
		}
	}

	/* Return the inside point of this interval: */
	if (outward){
		return xy_t(start.x + search_x * zl, start.y + search_y * zl);
	}
	return xy_t(start.x + search_x * zr, start.y + search_y * zr);

}


static path_xy_t refine_bounding_polygon(const AugmentedProj& proj,
                                         path_xy_t& boundary, double xmin,
                                         double xmax, double ymin, double ymax,
                                         double atol, double bisection_offset,
                                         double minimum_node_distance)
{
	if (boundary.size() < 2)
		return boundary;

	/* Connect the last point: */
	boundary.push_back(boundary[0]);


	auto on_boundary = [=](const xy_t& xy) -> bool {
		return xy.x == xmin || xy.x == xmax || xy.y == ymin || xy.y == ymax;
	};

	xy_t last = boundary[0];
	/* Start the actual polygon path from the first boundary node: */
	path_xy_t poly;
	poly.push_back(last);
	/* A stack that will hold the next node: */
	std::stack<xy_t> next;
	const double min_distance2 = minimum_node_distance * minimum_node_distance;
	const double bisection_offset2 = bisection_offset * bisection_offset;
	for (auto it = boundary.cbegin()+1; it != boundary.cend(); ++it){
		/* Iterative refinement: */
		next.push(*it);
		while (!next.empty()){
			/* The central point: */
			xy_t center(0.5*(next.top().x + last.x),
			            0.5*(next.top().y + last.y));

			/* If all three are on the boundary, need not consider further: */
			if (on_boundary(last) && on_boundary(next.top())
			    && on_boundary(center))
			{
				last = next.top();
				poly.push_back(last);
				next.pop();
				continue;
			}

			/* The direction of the segment: */
			const double dx = next.top().x - last.x;
			const double dy = next.top().y - last.y;

			/* Find a point just inside the boundary, within `atol` of it: */
			xy_t shifted(find_boundary(center, dy, -dx, proj, atol));

			/* Make sure that the point is inside the bounding box.
			 * If we have to force the point into the bounding box, add it.
			 * Otherwise apply the `bisection_offset` criterion also known
			 * from "paths.hpp" (`crop_and_refine`). */
			if (shifted.x > xmax || shifted.x < xmin || shifted.y > ymax
			    || shifted.y < ymin)
			{
				shifted = boundary_intersection(shifted, center, xmin, xmax,
				                                ymin, ymax);

				/* Add point: */
				if (last.distance_squared(shifted) <= min_distance2){
					/* Add point to polygon and continue. */
					poly.push_back(shifted);
					last = shifted;
					next.pop();
				} else {
					/* Add point to queue. */
					next.push(shifted);
				}
			} else {
				/* Apply the `bisection_offset` and minimum distance
				 * criteria. */
				if (center.distance_squared(shifted) <= bisection_offset2
				    || shifted.distance_squared(last) <= min_distance2)
				{
					/* Add point to polygon and continue. */
					poly.push_back(shifted);
					last = shifted;
					next.pop();
				} else {
					/* Add point to queue. */
					next.push(shifted);
				}
			}
		}
	}

	/* Disconnect the last point (undo what we added in the beginning): */
	if (poly.empty())
		throw std::runtime_error("Empty poly where it should not be.");
	poly.resize(poly.size()-1);


	return poly;
}


void flottekarte::bounding_polygon(const AugmentedProj&  proj, double xmin,
                                   double xmax, double ymin, double ymax,
                                   double atol, double bisection_offset,
                                   double minimum_node_distance,
                                   path_xy_t& poly)
{
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
			bool bij = inside(xy_ij, proj);
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
	if (boundary.empty()){

		/* Need to find the boundary. */
		boundary.push_back(find_boundary(xy_t(xmin, 0.5*(ymin+ymax)),
		                                 -1.0, 0.0, proj, atol));
		boundary.push_back(find_boundary(xy_t(0.5*(xmax+xmin), ymin),
		                                 0.0, -1.0, proj, atol));
		boundary.push_back(find_boundary(xy_t(xmax, 0.5*(ymin+ymax)),
		                                 1.0, 0.0, proj, atol));
		boundary.push_back(find_boundary(xy_t(0.5*(xmax+xmin), ymax),
		                                 0.0, 1.0, proj, atol));
	}
	path_xy_t refined(refine_bounding_polygon(proj, boundary, xmin, xmax,
	                                          ymin, ymax, atol,
	                                          bisection_offset,
	                                          minimum_node_distance));
	poly.swap(refined);

}
