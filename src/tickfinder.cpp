
#include <../include/tickfinder.hpp>
#include <../include/invert.hpp>

#include <functional>
#include <set>
#include <cmath>

#include <boost/math/tools/roots.hpp>
using boost::math::tools::eps_tolerance;
using boost::math::tools::toms748_solve;
using projplot::compute_ticks;
using projplot::gradient_descent_inverse_project;
using projplot::axis_t;
using projplot::tick_t;
using projplot::geo_t;
using projplot::geo_degrees_t;
using projplot::ProjWrapper;
using projplot::GriddedInverter;

#include <iostream>

struct integer_level {
	double x;
	long level;

	integer_level(double x, long level) : x(x), level(level)
	{};
};


static std::vector<integer_level>
compute_integer_levels(std::function<double(double)> fun, double x0, double x1)
{
	/* Root finding configuration: */
	const eps_tolerance<double> tol(8);

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
				return fun(x) - yt;
			};
			std::uintmax_t max_iter(40);
			try {
				std::pair<double,double> xtlr
				   = toms748_solve(rootfun, xl, xr, tol, max_iter);

				levels.emplace_back(0.5*(xtlr.first + xtlr.second), yt);
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
projplot::compute_ticks(const ProjWrapper& proj, const GriddedInverter& ginv,
                        axis_t ax, tick_t tick, double xmin, double xmax,
                        double ymin, double ymax, double tick_spacing)
{
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

	auto fun = [&](double z) -> double {
		/* Get the current location on the map (boundary): */
		xy_t xy(gen_xy(z));

		/* Invert: */
		geo_t lola0(ginv(xy));
		geo_t lola(gradient_descent_inverse_project(proj, xy, lola0.lambda,
		                                            lola0.phi));

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
		geo_t lola0(ginv(xy));
		if (tick == TICK_LON)
			lola0.lambda = il.level * tick_spacing;
		else
			lola0.phi = il.level * tick_spacing;
		geo_t lola(gradient_descent_inverse_project(proj, xy,
		                                            lola0.lambda,
		                                            lola0.phi));
		if (tick == TICK_LON)
			res.push_back({il.level * tick_spacing, rad2deg(lola.phi)});
		else
			res.push_back({rad2deg(lola.lambda), il.level * tick_spacing});

	}

	return res;
}
