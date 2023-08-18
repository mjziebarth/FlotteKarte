/*
 * Utility for finding ticks along an axis.
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
#include <../include/tickfinder.hpp>
#include <../include/invert.hpp>
#include <../include/gradient.hpp>

#include <functional>
#include <set>
#include <cmath>

#include <boost/math/tools/roots.hpp>
using boost::math::tools::eps_tolerance;
using boost::math::tools::toms748_solve;
using flottekarte::compute_ticks;
using flottekarte::gradient_descent_inverse_project;
using flottekarte::axis_t;
using flottekarte::tick_t;
using flottekarte::xy_t;
using flottekarte::geo_t;
using flottekarte::geo_degrees_t;
using flottekarte::segment_tick_t;
using flottekarte::ProjWrapper;
using flottekarte::GriddedInverter;
using flottekarte::Gradient;
using flottekarte::AugmentedProj;

#include <iostream>

segment_tick_t::segment_tick_t(size_t segment, double lon, double lat)
   : segment(segment), tick(lon, lat)
{
}

struct integer_level {
	double x;
	long level;

	integer_level(double x, long level) : x(x), level(level)
	{};
};

/* Root finding configuration: */
static const eps_tolerance<double> tol(26);

static std::vector<integer_level>
compute_integer_levels(std::function<double(double)> fun, double x0, double x1)
{
	/* Very simple grid search to identify all the rounded integers: */
	std::vector<integer_level> levels;
	const double dx = (x1 - x0) / 1000.0;
	std::vector<int> transitions;
	double xl = x0;
	double yl = fun(xl);
	double xr, yr;
	long yintl = std::floor(yl);
	for (int i=1; i<1000; ++i){
		/* A step forward: */
		xr = x0 + i*dx;
		yr = fun(xr);
		long yintr = std::floor(yr);
		if (yintl != yintr){
			/* Found a transition! */
			long yt = std::max(yintl,yintr);
			auto rootfun = [&](double x) -> double {
				return fun(x) - static_cast<double>(yt);
			};
			std::uintmax_t max_iter(100);
			try {
				std::pair<double,double> xtlr
				   = toms748_solve(rootfun, xl, xr, tol, max_iter);

				/* Check the integer level quality: */
				double zprop = 0.5 * (xtlr.first + xtlr.second);
				if (std::abs(rootfun(zprop)) < 1e-5){
					levels.emplace_back(0.5*(xtlr.first + xtlr.second), yt);
				}
			} catch (...) {
				/* Something went wrong. Better skip the tick than to
				 * panic Python. */
			}


			/* Proceed with the new level: */
			xl = xr;
			yl = yr;
			yintl = yintr;
		}
	}

	return levels;
}

static inline xy_t gen_xy(double z, double x0, double y0, double x1, double y1)
{
	return xy_t(z * x0 + (1.0 - z) * x1,   z * y0 + (1.0 - z) * y1);
}

static inline double cross_product_norm(double dx0, double dy0, double dx1,
                                        double dy1)
{
	/* Compute the norm of the cross product of two planar vectors. */
	return std::abs(dx0 * dy1 - dy0 * dx1);
}


std::vector<segment_tick_t>
flottekarte::compute_ticks(const AugmentedProj& proj,
                           tick_t tick, const path_xy_t& path,
                           double tick_spacing)
{
	/* No ticks: */
	if (tick == TICK_NONE || path.size() < 2)
		return std::vector<segment_tick_t>();

	/* Compute coordinate extremes: */
	double xmin = std::numeric_limits<double>::infinity();
	double xmax = -std::numeric_limits<double>::infinity();
	double ymin = std::numeric_limits<double>::infinity();
	double ymax = -std::numeric_limits<double>::infinity();
	for (const xy_t& xy : path){
		xmin = std::min(xmin, xy.x);
		xmax = std::max(xmax, xy.x);
		ymin = std::min(ymin, xy.y);
		ymax = std::max(ymax, xy.y);
	}

	std::vector<std::vector<segment_tick_t>> ticks(path.size());
	#pragma omp parallel for
	for (size_t i=0; i < path.size(); ++i){
		const size_t j = (i+1) % path.size();
		const double x0 = path[i].x, y0 = path[i].y,
		             x1 = path[j].x, y1 = path[j].y;
		const double segment_length = path[i].distance(path[j]);

		/* The function referring a point along the line segment to
		 * a coordinate: */
		auto fun = [&](double z) -> double {
			/* Get the current location on the map (boundary): */
			xy_t xy(gen_xy(z,x0,y0,x1,y1));

			/* Invert: */
			geo_t lola(proj.inverse(xy));

			/* Return the relevant coordinate scaled to spacing: */
			if (tick == TICK_LON){
				return rad2deg(lola.lambda) / tick_spacing;
			} else {
				return rad2deg(lola.phi) / tick_spacing;
			}
		};

		/* Compute: */
		std::vector<integer_level> int_levels
		   = compute_integer_levels(fun, 0.0, 1.0);

		/* Early exit if no levels: */
		if (int_levels.empty())
			continue;

	// /* Check whether the 180°=-180° meridian can be part of the
	//  * tick set: */
	// if (tick == TICK_LON){
	// 	for (size_t i=0; i<int_levels.size()-1; ++i){
	// 		if (((int_levels[i].level <= 0) != (int_levels[i+1].level <= 0))
	// 		    && (std::abs(int_levels[i].level - int_levels[i+1].level) > 2))
	// 		{
	// 			/* We cross the 180° meridian. Try to add it to the tick set: */
	// 			auto fun2 = [&](double phi) -> double {
	// 				if (ax == AX_BOT)
	// 					return proj.project(PI, phi).y - ymin;
	// 				else if (ax == AX_TOP)
	// 					return proj.project(PI, phi).y - ymax;
	// 				else if (ax == AX_LEFT)
	// 					return proj.project(PI, phi).x - xmin;
	// 				else if (ax == AX_RIGHT)
	// 					return proj.project(PI, phi).x - xmax;
	// 			};

	// 			/* First estimate of phi and finding a bracketing interval: */
	// 			const double z0 = 0.5 * (int_levels[i].x + int_levels[i+1].x);
	// 			const double phi0 = invert(gen_xy(z0)).phi;
	// 			double dp = 5e-6;
	// 			double phi_min, phi_max;
	// 			bool success = false;
	// 			for (int i=0; i<20; ++i){
	// 				dp *= 2.0;
	// 				phi_min = std::max(phi0 - dp, -PI/2);
	// 				phi_max = std::min(phi0 + dp, PI/2);
	// 				if ((fun2(phi_min) <= 0) != (fun2(phi_max) <= 0)){
	// 					success = true;
	// 					break;
	// 				}
	// 			}

	// 			if (success){
	// 				/* Use TOMS 748 to find the correct phi: */
	// 				std::uintmax_t max_iter = 40;
	// 				std::pair<double,double> phi180b
	// 				   = toms748_solve(fun2, phi_min, phi_max, tol, max_iter);
	// 				double phi180 = 0.5 * (phi180b.first + phi180b.second);

	// 				/* Add to ticks: */
	// 				int_levels.push_back({proj.project(PI, phi180).x,
	// 				                      std::round(180.0 / tick_spacing)});
	// 			}

	// 			break;
	// 		}
	// 	}
	// }

		/* Here we sort out some possible duplicates.
		 * We collect all of the ticks within `int_levels`
		 * into vectors of the same long `level`.
		 * Within these vector bins, we collect only 'unique' ticks,
		 * that is, ticks which are further than 1e-2*(zmax-zmin)
		 * apart in map coordinates.
		 */
		const double min_spacing = 1e-2*std::max(xmax-xmin, ymax-ymin);
		std::map<long,std::vector<double>> multilevels;
		for (auto it = int_levels.begin(); it != int_levels.end(); ++it) {
			auto find = multilevels.find(it->level);
			if (find != multilevels.cend()){
				/* There is already at least one tick of that level.
				 * Check all existing links: */
				bool unique = true;
				for (double x : find->second){
					if (std::abs(x - it->x) * segment_length < min_spacing){
						/* Another very close tick exists. Not unique tick,
						 * keep the other! */
						unique = false;
						break;
					}
				}
				if (unique){
					/* Keep in list. */
					find->second.push_back(it->x);
				}
			} else {
				/* First or unique one: */
				auto it2 = multilevels.emplace(it->level, 0).first;
				it2->second.push_back(it->x);
			}
		}
		int_levels.clear();
		for (auto ticks : multilevels) {
			for (double x : ticks.second) {
				int_levels.emplace_back(x, ticks.first);
			}
		}


		//std::vector<geo_degrees_t> res;
		//res.reserve(int_levels.size());
		for (auto il : int_levels){
			xy_t xy(gen_xy(il.x, x0, y0, x1, y1));
			geo_t lola(proj.inverse(xy));

			/* Sanity check: Make sure that a full projection-inversion loop
			 * can be performed. */
			double xy_norm = std::sqrt(xy.x*xy.x + xy.y * xy.y);
			if (proj.project(lola).distance(xy) > 1e-5*xy_norm)
				continue;

			/* Make sure that the tick is not too flat relative to
			 * the segment: */
			constexpr double sin_angle_min = std::sin(10.0 * PI/180.0);
			Gradient<FORWARD_5POINT> gr(proj.wrapper(), lola);

			const double dx0 = (x1 - x0) / segment_length;
			const double dy0 = (y1 - y0) / segment_length;
			if (tick == TICK_LON){
				const double norm1 = std::sqrt(  gr.gx_north() * gr.gx_north()
				                               + gr.gy_north() * gr.gy_north());
				const double dx1 = gr.gx_north() / norm1;
				const double dy1 = gr.gy_north() / norm1;
				const double sin_angle \
				   = cross_product_norm(dx0, dy0, dx1, dy1);
				if (sin_angle >= sin_angle_min)
					ticks[i].emplace_back(i, il.level * tick_spacing,
					                      rad2deg(lola.phi));
			} else {
				const double norm1 = std::sqrt(  gr.gx_east() * gr.gx_east()
				                               + gr.gy_east() * gr.gy_east());
				const double dx1 = gr.gx_east() / norm1;
				const double dy1 = gr.gy_east() / norm1;
				const double sin_angle \
				   = cross_product_norm(dx0, dy0, dx1, dy1);
				if (sin_angle >= sin_angle_min)
					ticks[i].emplace_back(i, rad2deg(lola.lambda),
					                      il.level * tick_spacing);
			}
		}
	}

	size_t Nticks = 0;
	for (auto& vec : ticks)
		Nticks += vec.size();
	std::vector<segment_tick_t> ticks_serial;
	ticks_serial.reserve(Nticks);
	for (auto& vec : ticks){
		ticks_serial.insert(ticks_serial.cend(), vec.cbegin(), vec.cend());
		vec.clear();
	}

	return ticks_serial;
}
