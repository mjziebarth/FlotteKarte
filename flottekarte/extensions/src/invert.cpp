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

#include <../include/invert.hpp>
#include <../include/linalg.hpp>
#include <iostream>

#include <boost/math/tools/roots.hpp>
using boost::math::tools::eps_tolerance;
using boost::math::tools::bracket_and_solve_root;

using flottekarte::geo_t;
using flottekarte::ProjWrapper;
using flottekarte::xy_t;

flottekarte::geo_t
flottekarte::gradient_descent_inverse_project(const ProjWrapper& projection,
                                              const xy_t& xy, double lambda0,
                                              double phi0)
{
	/* Use root finding to compute inverse. */

	/* Problem definition: */
	typedef linalg_t<2,double> lina;
	typedef typename lina::column_vectord_t vd_t;
	vd_t lola;
	lola[0] = lambda0;
	lola[1] = phi0;

	/* Wrapping lambda and phi continuously: */
	bool have_warned = false;
	auto compute_lambda_phi = [&](const vd_t& la_ph) -> geo_t {
		if (la_ph[1] > 1e3*PI || la_ph[1] < -1e3*PI){
			// if (!have_warned){
			// 	std::cerr << "Warning: |phi| > 1e3 * pi!\n";
			// 	std::cerr << "  lambda =" << la_ph[0] << "\n";
			// 	std::cerr << "     phi =" << la_ph[1] << "\n";
			// 	have_warned = true;
			// }
		}
		geo_t geo;
		const int winding_number = (la_ph[1] >= 0.0) ?
		                              std::floor(std::abs(la_ph[1]+PI/2) / PI)
		                           : std::floor(std::abs(la_ph[1]-PI/2) / PI);
		const bool even = (winding_number % 2) == 0;
		geo.lambda = (even) ? modulo(la_ph[0], 2*PI) : modulo(-la_ph[0],2*PI);
		geo.phi = modulo(la_ph[1]+PI/2, PI) - PI/2;
		if (!even)
			geo.phi = -geo.phi;
		if (geo.phi > PI/2)
			std::cerr << "geo.phi > PI/2\n" << std::flush;
		if (geo.phi < -PI/2)
			std::cerr << "geo.phi < -PI/2\n" << std::flush;
		return geo;
	};

	auto cost = [&](const vd_t& lola) -> double {
		const geo_t geo(compute_lambda_phi(lola));
		const xy_t xy_proj(projection.project(geo.lambda, geo.phi));
		const double dx = (xy_proj.x - xy.x);
		const double dy = (xy_proj.y - xy.y);
		return dx * dx + dy * dy;
	};

	double cost_i = cost(lola);
	auto gradient = [&](const vd_t& x) -> vd_t {
		/* Numerical gradient. */
		constexpr double delta = 1e-5;
		const vd_t dx({delta, 0.0});
		const vd_t dy({0.0, delta});
		vd_t g;
		/* Third-order accurate forward finite difference: */
		constexpr double f3 = 1.0/3.0, f2 = -1.5, f1 = 3.0, f0 = 11.0/6.0;
		g[0] = 1.0/delta * (f3*cost(x+3*dx) + f2*cost(x+2*dx)
		                    + f1*cost(x+dx) - f0 * cost_i);
		g[1] = 1.0/delta * (f3*cost(x+3*dy) + f2*cost(x+2*dy)
		                    + f1*cost(x+dy) - f0 * cost_i);
		return g;
	};


	constexpr int MAX_STEPS = 1000;
	const eps_tolerance<double> tol(12);
	double step = 0.1;
	int i, m=0;
	for (i=0; i<MAX_STEPS; ++i){
		if (cost_i < 1e-30){
			break;
		}

		/* Compute gradient: */
		vd_t g(gradient(lola));

		/* Every ever-so-many iterations, check whether we are in
		 * a valley of death close to the poles.
		 * There, we might have a super slim valley: gradients in
		 * latitude direction are very steep but the longitude
		 * (which remains indetermined) direction is very shallow. */
		if (m >= 20){
			auto cost_lon = [&](double lambda)->double {
				const geo_t geo(compute_lambda_phi(vd_t({lambda,
						                                 lola[1]})));
				const xy_t xy_proj(projection.project(geo.lambda,
						                              geo.phi));
				const double dx = (xy_proj.x - xy.x);
				const double dy = (xy_proj.y - xy.y);
				return dx + dy;
			};

			/* Bring back to coordinate plane: */
			geo_t geo(compute_lambda_phi(lola));
			lola[0] = geo.lambda;
			lola[1] = geo.phi;
			bool is_rising = cost_lon(lola[0]+1e-5)
			                    > cost_lon(lola[0]);
			std::uintmax_t iterations = 40;
			try {
				std::pair<double,double> lr
				   = bracket_and_solve_root(cost_lon, lola[0], 1.3,
				                            is_rising, tol, iterations);
				vd_t lola_propose({0.5 * (lr.first + lr.second),
				                   lola[1]});
				double cost_propose = cost(lola_propose);
				if (cost_propose < cost_i){
					/* Accept. */
					lola = lola_propose;
					cost_i = cost_propose;
				}
			} catch (...) {
				/* Failed the root finding, reset counter and continue
				 * as if nothing ever happened (except the parameter
				 * cropping). */
			}

			/* Reset the counter: */
			m = 0;

			/* Do not continue as usual: */
			cost_i = cost(lola);
			continue;
		} else {
			++m;
		}

		/* Make sure that we choose a step size that reduces the cost: */
		double cost_n0 = cost(lola - step*g);
		while (cost_n0 > cost_i && lina::norm(step*g) > 1e-16){
			step *= 0.125;
			cost_n0 = cost(lola - step*g);
		}

		if (lina::norm(step*g) <= 1e-16){
			break;
		}

		/* Compare costs for different variations of step size: */
		constexpr double ADJ = 1.5;
		if ((i % 2) == 0){
			const double cost_n1 = cost(lola - ADJ*step*g);
			if (cost_n1 < cost_n0){
				step *= ADJ;
				cost_i = cost_n1;
			} else {
				cost_i = cost_n0;
			}
		} else {
			const double cost_n1 = cost(lola - 1.0/ADJ*step*g);
			if (cost_n1 < cost_n0){
				step *= 1.0/ADJ;
				cost_i = cost_n1;
			} else {
				cost_i = cost_n0;
			}
		}

		/* Proceed in lon & lat: */
		lola -= step*g;
	}

	/* Polish the coordinate results: */
	geo_t lp(compute_lambda_phi(lola));


	if (lp.lambda > PI)
		lp.lambda -= 2*PI;

	return lp;
}
