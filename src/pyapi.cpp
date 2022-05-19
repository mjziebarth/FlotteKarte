/*
 *
 */

#include <../include/pyapi.hpp>
#include <../include/invert.hpp>
#include <../include/projwrapper.hpp>
#include <../include/griddedinverter.hpp>
#include <../include/gradient.hpp>
#include <../include/tickfinder.hpp>
#include <../include/grid.hpp>
#include <iostream>

using projplot::xy_t;
using projplot::geo_t;
using projplot::axis_t;
using projplot::tick_t;
using projplot::rad2deg;
using projplot::deg2rad;
using projplot::path_xy_t;
using projplot::ProjError;
using projplot::ProjWrapper;
using projplot::geo_degrees_t;
using projplot::GriddedInverter;
using projplot::generate_grid_lines;
using projplot::gradient_descent_inverse_project;
using projplot::compute_ticks;
using projplot::Gradient;
using projplot::FORWARD_5POINT;

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
			try {
				/* First try to use PROJ inverse of projection: */
				geo_t lola(proj.inverse(xy));
				lon_lat[2*i] = rad2deg(lola.lambda);
				lon_lat[2*i+1] = rad2deg(lola.phi);
			} catch (const ProjError& e) {
				/* If this fails, invert with gradient descent: */
				geo_t start(ginv(xy));
				geo_t lola(gradient_descent_inverse_project(proj, xy,
				                                            start.lambda,
				                                            start.phi));
				lon_lat[2*i] = rad2deg(lola.lambda);
				lon_lat[2*i+1] = rad2deg(lola.phi);
			}
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

int gradients_east_north(const char* proj_str, unsigned long Npoints,
                          const double* lon, const double* lat,
                          double* gradient_east, double* gradient_north,
                          double stencil_delta)
{
	try {
		/* Initialize the projection: */
		ProjWrapper proj(proj_str);

		/* Parallel evaluation of the gradients: */
		#pragma omp parallel for
		for (size_t i=0; i<Npoints; ++i){
			/* Compute coordinate gradients in east and north: */
			geo_t lola({deg2rad(lon[i]), deg2rad(lat[i])});
			Gradient<FORWARD_5POINT> gradient(proj, lola, stencil_delta);

			/* Save the values: */
			gradient_east[2*i]    = gradient.gx_east();
			gradient_east[2*i+1]  = gradient.gy_east();
			gradient_north[2*i]   = gradient.gx_north();
			gradient_north[2*i+1] = gradient.gy_north();
		}

		/* Success: */
		return 0;
	} catch (...) {
		return 1;
	}
}

int scale_factors(const char* proj_str, unsigned long Npoints,
                  const double* lon, const double* lat,
                  double* kh, double stencil_delta)
{
	try {
		/* Initialize the projection: */
		ProjWrapper proj(proj_str);

		/* Parallel evaluation of the gradients: */
		#pragma omp parallel for
		for (size_t i=0; i<Npoints; ++i){
			/* Compute coordinate gradients in east and north: */
			geo_t lola({deg2rad(lon[i]), deg2rad(lat[i])});
			Gradient<FORWARD_5POINT> g(proj, lola, stencil_delta);

			/* Save the values: */
			const double k = std::sqrt(  g.gx_east() * g.gx_east()
			                           + g.gy_east() * g.gy_east());
			kh[2*i]   = k;
			const double h = std::sqrt(  g.gx_north() * g.gx_north()
			                           + g.gy_north() * g.gy_north());
			kh[2*i+1] = h;
		}

		/* Success: */
		return 0;
	} catch (...) {
		return 1;
	}
}

int compute_axes_ticks(const char* proj_str, double xmin, double xmax,
                       double ymin, double ymax, double tick_spacing_degree,
                       unsigned int max_ticks_per_axis,
                       unsigned char which_bot, unsigned char which_top,
                       unsigned char which_left, unsigned char which_right,
                       double* bot_ticks, double* top_ticks,
                       double* left_ticks, double* right_ticks,
                       unsigned int* Nticks)
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
		std::cout << "tick creation loop...\n" << std::flush;
		#endif

		/* Prepare everything for loop execution: */
		std::array<axis_t,4> axes({projplot::AX_BOT,
		                           projplot::AX_TOP,
		                           projplot::AX_LEFT,
		                           projplot::AX_RIGHT});
		std::array<tick_t,4> tick({static_cast<tick_t>(which_bot),
		                           static_cast<tick_t>(which_top),
		                           static_cast<tick_t>(which_left),
		                           static_cast<tick_t>(which_right)});
		std::array<double*,4> tick_out({bot_ticks, top_ticks, left_ticks,
		                                right_ticks});
		std::array<std::vector<geo_t>,4> ticks;

		const size_t Nmax = static_cast<size_t>(max_ticks_per_axis);
		#pragma omp parallel for
		for (int i=0; i<4; ++i){
			std::vector<geo_degrees_t> ticks
			   = compute_ticks(proj, ginv, axes[i], tick[i], xmin, xmax,
			                   ymin, ymax, tick_spacing_degree);
			size_t Ni = std::min<size_t>(ticks.size(), Nmax);
			for (size_t j=0; j<Ni; ++j){
				*(tick_out[i] + 2*j) = ticks[j].lon;
				*(tick_out[i] + 2*j + 1) = ticks[j].lat;
			}
			Nticks[i] = static_cast<unsigned int>(Ni);
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

struct grid_lines_t {
	std::vector<path_xy_t> paths;
};

/*
 * This is the first of a two-part function. Computes grid lines according to
 * settings (projection, x- & ylim, tick spacing, bisection offset,
 * minimum distance between path-adjacent nodes) and returns 1) the length
 * of the resulting path (Npath) and 2) a pointer to the structure holding
 * all the required information (struct_ptr).
 * It has to be followed, after allocating two suitable numpy arrays,
 * by a call to save_grid_lines, which transfers the path contained in
 * struct_ptr and frees the allocated space. If this call is omitted, a memory
 * leak will occur. If this call a second time, unallocated memory will be
 * accessed.
 */
int compute_grid_lines(const char* proj_str, double xmin, double xmax,
                       double ymin, double ymax, int tick_spacing_degree,
                       double bisection_offset, double minimum_node_distance,
                       double max_lat, double cut_at_angle_degrees,
                       void** struct_ptr, size_t* Npath)
{
	/* Empty initialization: */
	if (!struct_ptr){
		std::cerr << "No pointer to write struct given.\n";
		return 3;
	}
	if (!Npath){
		std::cerr << "No pointer to write Npath given.\n";
		return 4;
	}
	*struct_ptr = nullptr;
	*Npath = 0;

	try {
		/* Create the projection wrapper: */
		ProjWrapper proj(proj_str);

		/* Generate the grid lines structure: */
		grid_lines_t* glines = new grid_lines_t();

		/* Save it to the output pointer: */
		*struct_ptr = static_cast<void*>(glines);

		/* Compute the grid lines: */
		try {
			glines->paths = generate_grid_lines(proj, xmin, xmax, ymin, ymax,
			                                    tick_spacing_degree,
			                                    bisection_offset,
			                                    minimum_node_distance, max_lat,
			                                    cut_at_angle_degrees);

		} catch (...) {
			/* Clean up: */
			delete glines;
			struct_ptr = nullptr;

			return 1;
		}

		/* Compute the number of grid lines: */
		size_t npath = 0;
		for (const path_xy_t& path : glines->paths){
			npath += path.size();
		}
		*Npath = npath;

		/* Success. */
		return 0;

	} catch (const ProjError& err){
		/* Could not create the projection object: */
		std::cerr << "ProjError: " << err.what() << "\n";
		return 2;
	}
}

constexpr uint8_t MPL_MOVETO = 1;
constexpr uint8_t MPL_LINETO = 2;

/*
 * Second part of a two-part function. Call exactly one time after
 * compute_grid_lines has successfully completed.
 */
int save_grid_lines(const void* struct_ptr, double* vertices, uint8_t* codes)
{
	/* Cast the struct: */
	if (struct_ptr == nullptr)
		return 1;
	const grid_lines_t& glines = *static_cast<const grid_lines_t*>(struct_ptr);

	/* Fill the vertices and codes arrays: */
	for (const path_xy_t& path : glines.paths){
		for (size_t i=0; i<path.size(); ++i){
			/* Set the code: */
			if (i == 0){
				*codes = MPL_MOVETO;
			} else {
				*codes = MPL_LINETO;
			}

			/* Set the vertices: */
			*vertices = path[i].x;
			++vertices;
			*vertices = path[i].y;

			/* Advance: */
			++codes;
			++vertices;
		}
	}

	/* Success. */
	return 0;
}

int clean_grid_lines_struct(void* struct_ptr)
{
	if (struct_ptr == nullptr)
		return 1;

	/* Delete the grid_lines_t struct, calling all necessary destructors: */
	grid_lines_t* glines = static_cast<grid_lines_t*>(struct_ptr);
	delete glines;

	/* Success. */
	return 0;
}
