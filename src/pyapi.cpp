/*
 * API used for interfacing with Python.
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

#include <../include/pyapi.hpp>
#include <../include/invert.hpp>
#include <../include/projwrapper.hpp>
#include <../include/griddedinverter.hpp>
#include <../include/gradient.hpp>
#include <../include/tickfinder.hpp>
#include <../include/grid.hpp>
#include <iostream>

using flottekarte::xy_t;
using flottekarte::geo_t;
using flottekarte::cut_t;
using flottekarte::axis_t;
using flottekarte::tick_t;
using flottekarte::rad2deg;
using flottekarte::deg2rad;
using flottekarte::path_xy_t;
using flottekarte::ProjError;
using flottekarte::geo_grid_t;
using flottekarte::ProjWrapper;
using flottekarte::geo_degrees_t;
using flottekarte::GriddedInverter;
using flottekarte::generate_grid_lines;
using flottekarte::gradient_descent_inverse_project;
using flottekarte::compute_ticks;
using flottekarte::Gradient;
using flottekarte::FORWARD_5POINT;

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
		std::array<axis_t,4> axes({flottekarte::AX_BOT,
		                           flottekarte::AX_TOP,
		                           flottekarte::AX_LEFT,
		                           flottekarte::AX_RIGHT});
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
	std::vector<cut_t> cuts;
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
                       double max_lat, void** struct_ptr, size_t* Npath,
                       size_t* Ncut)
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
	if (!Ncut){
		std::cerr << "No pointer to write Ncut given.\n";
		return 5;
	}
	*struct_ptr = nullptr;
	*Npath = 0;
	*Ncut = 0;

	try {
		/* Create the projection wrapper: */
		ProjWrapper proj(proj_str);

		/* Generate the grid lines structure: */
		grid_lines_t* glines = new grid_lines_t();

		/* Save it to the output pointer: */
		*struct_ptr = static_cast<void*>(glines);

		/* Compute the grid lines: */
		try {
			geo_grid_t grid(generate_grid_lines(proj, xmin, xmax, ymin, ymax,
			                                    tick_spacing_degree,
			                                    bisection_offset,
			                                    minimum_node_distance,
			                                    max_lat));
			glines->paths.swap(grid.paths);
			glines->cuts.swap(grid.cuts);

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

		/* Number of cuts: */
		*Ncut = glines->cuts.size();

		/* Success. */
		return 0;

	} catch (const ProjError& err){
		/* Could not create the projection object: */
		std::cerr << "ProjError: " << err.what() << "\n";
		return 2;
	} catch (const std::bad_alloc& err){
		/* Possibly could not allocate the grid lines object. */
		std::cerr << "Memory allocation error.\n";
		return 6;
	}
}

constexpr uint8_t MPL_MOVETO = 1;
constexpr uint8_t MPL_LINETO = 2;

/*
 * Second part of a two-part function. Call exactly one time after
 * compute_grid_lines has successfully completed.
 */
int save_grid_lines(const void* struct_ptr, double* vertices, uint8_t* codes,
                    double* cut_vertices, uint8_t* cut_axes, double* cut_coords)
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

	/* Fill the cut points: */
	for (const cut_t& cut : glines.cuts){
		*cut_vertices = cut.point.x;
		++cut_vertices;
		*cut_vertices = cut.point.y;
		++cut_vertices;
	}
	for (const cut_t& cut : glines.cuts){
		*cut_axes = static_cast<uint8_t>(cut.tick_type);
		++cut_axes;
	}
	for (const cut_t& cut : glines.cuts){
		*cut_coords = cut.coordinate;
		++cut_coords;
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
