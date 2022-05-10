
#include <../include/invert.hpp>
#include <../include/linalg.hpp>
#include <iostream>

using projplot::geo_t;
using projplot::ProjWrapper;
using projplot::xy_t;

projplot::geo_t
projplot::gradient_descent_inverse_project(const ProjWrapper& projection,
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
		return geo;
	};

	auto cost = [&](const vd_t& lola) -> double {
		const geo_t geo(compute_lambda_phi(lola));
		const xy_t xy_proj(projection.project(geo.lambda, geo.phi));
		const double dx = (xy_proj.x - xy.x);
		const double dy = (xy_proj.y - xy.y);
		return dx * dx + dy * dy;
	};

	auto gradient = [&](const vd_t& x) -> vd_t {
		/* Numerical gradient. */
		constexpr double delta = 1e-5;
		const vd_t dx({delta, 0.0});
		const vd_t dy({0.0, delta});
		vd_t g;
		/* Fourth-order accurate symmetric difference: */
		g[0] = 1.0/(12*delta) * (-cost(x+2*dx) + 8*cost(x+dx)
		                         - 8*cost(x-dx) + cost(x-2*dx));
		g[1] = 1.0/(12*delta) * (-cost(x+2*dy) + 8*cost(x+dy)
		                         - 8*cost(x-dy) + cost(x-2*dy));
		return g;
	};


	double step = 0.5;
	double cost_i = cost(lola);
	int i;
	for (i=0; i<1000; ++i){
		if (cost_i < 1e-30){
			break;
		}

		/* Compute gradient: */
		vd_t g(gradient(lola));

		/* Make sure that we choose a step size that reduces the cost: */
		double cost_n0 = cost(lola - step*g);
		while (cost_n0 > cost_i && lina::norm(step*g) > 1e-16){
			step *= 0.5;
			cost_n0 = cost(lola - step*g);
		}

		if (lina::norm(step*g) <= 1e-16){
			break;
		}

		/* Compare costs for different variations of step size: */
		const double cost_n1 = cost(lola - 1.1*step*g);
		const double cost_n2 = cost(lola - 0.9*step*g);
		if (cost_n1 < cost_n0){
			if (cost_n1 < cost_n2){
				cost_i = cost_n1;
				step *= 1.1;
			} else {
				cost_i = cost_n2;
				step *= 0.9;
			}
		} else if (cost_n2 < cost_n0){
			cost_i = cost_n1;
			step *= 0.9;
		} else {
			cost_i = cost_n0;
		}

		/* Proceed in lon & lat: */
		lola -= step*g;
		//std::cout << "    step[" << i << "]: lola = (" << lola[0] << ","
		//          << lola[1] << ")\n";
	}

	/* Polish the coordinate results: */
	geo_t lp(compute_lambda_phi(lola));
	//std::cout << "got (" << lp.lambda << "," << lp.phi << ")\n" << std::flush;


	if (lp.lambda > PI)
		lp.lambda -= 2*PI;

	return lp;
}
