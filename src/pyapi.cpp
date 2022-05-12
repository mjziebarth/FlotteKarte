/*
 *
 */

#include <../include/pyapi.hpp>
#include <../include/invert.hpp>
#include <../include/projwrapper.hpp>
#include <../include/griddedinverter.hpp>
#include <iostream>

using projplot::xy_t;
using projplot::geo_t;
using projplot::rad2deg;
using projplot::deg2rad;
using projplot::ProjError;
using projplot::ProjWrapper;
using projplot::GriddedInverter;
using projplot::gradient_descent_inverse_project;

/*void project_data(const char* proj_str, unsigned long Npoints,
                  const double* lon, const double* lat,
                  double* xy)
{

}*/

int inverse_project_data_optimize(const char* proj_str, unsigned long Npoints,
                                   const double* x, const double* y,
                                   double* lon_lat)
{
	/* Initialize the projection: */
	#ifdef DEBUG
	std::cout << "initializing proj wrapper\n" << std::flush;
	#endif
	try {
		ProjWrapper proj(proj_str);
		GriddedInverter ginv(proj, 100, 50);

		/* Parallel projection: */
		#ifdef DEBUG
		std::cout << "start projection loop...\n" << std::flush;
		#endif

		#pragma omp parallel for
		for (size_t i=0; i<Npoints; ++i){
			xy_t xy({x[i], y[i]});
			geo_t start(ginv(xy));
			geo_t lola(gradient_descent_inverse_project(proj, xy, start.lambda,
			                                            start.phi));
			lon_lat[2*i] = rad2deg(lola.lambda);
			lon_lat[2*i+1] = rad2deg(lola.phi);
		}

		#ifdef DEBUG
		std::cout << "cleaning up!\n" << std::flush;
		#endif

		/* Success: */
		return 0;
	} catch (const ProjError& err) {
		std::cout << "error occurred.\n what: " << err.what() << "\n"
		          << std::flush;
		return 1;
	} catch (...) {
		#ifdef DEBUG
		std::cout << "error occurred.\n" << std::flush;
		#endif
		return 1;
	}
}
