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
using flottekarte::geo_t;
using flottekarte::geo_degrees_t;
using flottekarte::ProjWrapper;
using flottekarte::GriddedInverter;

#include <iostream>

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


std::vector<geo_degrees_t>
flottekarte::compute_ticks(const ProjWrapper& proj, const GriddedInverter& ginv,
                        axis_t ax, tick_t tick, double xmin, double xmax,
                        double ymin, double ymax, double tick_spacing)
{
	/* No ticks: */
	if (tick == TICK_NONE)
		return std::vector<geo_degrees_t>();

	auto gen_xy = [&](double z) -> xy_t {
		xy_t xy;
		if (ax == AX_BOT){
			xy.x = z;
			xy.y = ymin;
		} else if (ax == AX_TOP){
			xy.x = z;
			xy.y = ymax;
		} else if (ax == AX_LEFT){
			xy.x = xmin;
			xy.y = z;
		} else {
			xy.x = xmax;
			xy.y = z;
		}
		return xy;
	};

	auto invert = [&](const xy_t& xy) -> geo_t {
		geo_t lola;
		try {
			lola = proj.inverse(xy);
		} catch (const ProjError& e) {
			geo_t lola0(ginv(xy));
			lola = gradient_descent_inverse_project(proj, xy, lola0.lambda,
			                                        lola0.phi);
		}
		return lola;
	};

	auto fun = [&](double z) -> double {
		/* Get the current location on the map (boundary): */
		xy_t xy(gen_xy(z));

		/* Invert: */
		geo_t lola(invert(xy));

		/* Return the relevant coordinate scaled to spacing: */
		if (tick == TICK_LON){
			return rad2deg(lola.lambda) / tick_spacing;
		} else {
			return rad2deg(lola.phi) / tick_spacing;
		}
	};

	/* Limits: */
	double zmin,zmax;
	if (ax == AX_BOT || ax == AX_TOP){
		zmin = xmin;
		zmax = xmax;
	} else {
		zmin = ymin;
		zmax = ymax;
	}

	/* Compute: */
	std::vector<integer_level> int_levels
	   = compute_integer_levels(fun, zmin, zmax);

	auto ax_label = [](axis_t ax) -> const char* {
		if (ax == AX_BOT)
			return "bot";
		else if (ax == AX_TOP)
			return "top";
		else if (ax == AX_LEFT)
			return "left";
		else if (ax == AX_RIGHT)
			return "right";
	};

	/* Check whether the 180°=-180° meridian can be part of the
	 * tick set: */
	if (tick == TICK_LON){
		for (size_t i=0; i<int_levels.size()-1; ++i){
			if (((int_levels[i].level <= 0) != (int_levels[i+1].level <= 0))
			    && (std::abs(int_levels[i].level - int_levels[i+1].level) > 2))
			{
				/* We cross the 180° meridian. Try to add it to the tick set: */
				auto fun2 = [&](double phi) -> double {
					if (ax == AX_BOT)
						return proj.project(PI, phi).y - ymin;
					else if (ax == AX_TOP)
						return proj.project(PI, phi).y - ymax;
					else if (ax == AX_LEFT)
						return proj.project(PI, phi).x - xmin;
					else if (ax == AX_RIGHT)
						return proj.project(PI, phi).x - xmax;
				};

				/* First estimate of phi and finding a bracketing interval: */
				const double z0 = 0.5 * (int_levels[i].x + int_levels[i+1].x);
				const double phi0 = invert(gen_xy(z0)).phi;
				double dp = 5e-6;
				double phi_min, phi_max;
				bool success = false;
				for (int i=0; i<20; ++i){
					dp *= 2.0;
					phi_min = std::max(phi0 - dp, -PI/2);
					phi_max = std::min(phi0 + dp, PI/2);
					if ((fun2(phi_min) <= 0) != (fun2(phi_max) <= 0)){
						success = true;
						break;
					}
				}

				if (success){
					/* Use TOMS 748 to find the correct phi: */
					std::uintmax_t max_iter = 40;
					std::pair<double,double> phi180b
					   = toms748_solve(fun2, phi_min, phi_max, tol, max_iter);
					double phi180 = 0.5 * (phi180b.first + phi180b.second);

					/* Add to ticks: */
					int_levels.push_back({proj.project(PI, phi180).x,
					                      std::round(180.0 / tick_spacing)});
				}

				break;
			}
		}
	}

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
				if (std::abs(x - it->x) < min_spacing){
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


	std::vector<geo_degrees_t> res;
	res.reserve(int_levels.size());
	for (auto il : int_levels){
		xy_t xy(gen_xy(il.x));
		geo_t lola(invert(xy));

		/* Sanity check: Make sure that a full projection-inversion loop
		 * can be performed. */
		double xy_norm = std::sqrt(xy.x*xy.x + xy.y * xy.y);
		if (proj.project(lola).distance(xy) > 1e-5*xy_norm)
			continue;

		if (tick == TICK_LON)
			res.push_back({il.level * tick_spacing, rad2deg(lola.phi)});
		else
			res.push_back({rad2deg(lola.lambda), il.level * tick_spacing});

	}

	return res;
}
