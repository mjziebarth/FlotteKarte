/*
 *
 */

#include <../include/pyapi.hpp>
#include <../include/invert.hpp>
#include <../include/projwrapper.hpp>
#include <iostream>

using projplot::xy_t;
using projplot::geo_t;
using projplot::rad2deg;
using projplot::deg2rad;
using projplot::ProjWrapper;
using projplot::gradient_descent_inverse_project;

/*void project_data(const char* proj_str, unsigned long Npoints,
                  const double* lon, const double* lat,
                  double* xy)
{

}*/

void inverse_project_data_optimize(const char* proj_str, unsigned long Npoints,
                                   const double* x, const double* y,
                                   double lon0, double lat0,
                                   double* lon_lat)
{
	/* Initialize the projection: */
	std::cout << "  lon0 = " << lon0 << "\n";
	std::cout << "  lat0 = " << lat0 << "\n";
	std::cout << "initializing proj wrapper\n" << std::flush;
	try {
		ProjWrapper proj(proj_str);
		
		/* Parallel projection: */
		const double lambda0 = deg2rad(lon0);
		const double phi0 = deg2rad(lat0);
		std::cout << "start projection loop...\n" << std::flush;
		#pragma omp parallel for
		for (size_t i=0; i<Npoints; ++i){
			xy_t xy({x[i], y[i]});
			geo_t lola(gradient_descent_inverse_project(proj, xy, lambda0, phi0));
			lon_lat[2*i] = rad2deg(lola.lambda);
			lon_lat[2*i+1] = rad2deg(lola.phi);
		}
		std::cout << "cleaning up!\n" << std::flush;
	} catch (...) {
		std::cout << "error occurred.\n" << std::flush;
	}
}
