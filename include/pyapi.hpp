/*
 * API used for interfacing with Python.
 *
 * Authors: Malte J. Ziebarth (ziebarth@gfz-potsdam.de)
 *
 * Copyright (C) 2022 Deutsches GeoForschungsZentrum Potsdam,
 *                    Malte J. Ziebarth
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
 * */

#include <stdint.h>
#include <cstddef>

#ifndef FLOTTEKARTE_PYAPI_HPP
#define FLOTTEKARTE_PYAPI_HPP

extern "C" {

void project_data(const char* proj_str, unsigned long Npoints,
                  const double* lon, const double* lat,
                  double* xy);

void inverse_project_data(const char* proj_str, unsigned long Npoints,
                          const double* x, const double* y,
                          double* xy);

int inverse_project_data_optimize(const char* proj_str, unsigned long Npoints,
                                  const double* x, const double* y,
                                  double* xy);

int gradients_east_north(const char* proj_str, unsigned long Npoints,
                         const double* lon, const double* lat,
                         double* gradient_east, double* gradient_north,
                         double stencil_delta);

int scale_factors(const char* proj_str, unsigned long Npoints,
                  const double* lon, const double* lat,
                  double* kh, double stencil_delta);

int compute_axes_ticks(const char* proj_str, size_t Nseg,
                       const double* vertices, double tick_spacing_degree,
                       unsigned int max_ticks, unsigned int* segments,
                       double* tick_vertices, unsigned char* which_ticks,
                       unsigned int* Nticks);

/*
 * This is the first of a three-part function. Computes grid lines according to
 * settings (projection, x- & ylim, tick spacing, bisection offset,
 * minimum distance between path-adjacent nodes) and returns 1) the length
 * of the resulting path (Npath), 2) a pointer to the structure holding
 * all the required information (struct_ptr), and 3) the number of cuts (Ncut).
 * It should to be followed, after allocating three suitable numpy arrays,
 * by a call to save_grid_lines, which transfers the path contained in
 * struct_ptr and frees the allocated space.
 * It then has to be followed by a call to clean_grid_lines_struct. If this
 * last call is omitted, a memory leak will occur. If this is called a second
 * time, deallocated memory will be accessed.
 */
int compute_grid_lines(const char* proj_str, double xmin, double xmax,
                       double ymin, double ymax, int tick_spacing_degree,
                       double bisection_offset, double minimum_node_distance,
                       double max_lat, void** struct_ptr, size_t* Npath,
                       size_t* Ncut);

/*
 * Second part of a three-part function.
 */
int save_grid_lines(const void* struct_ptr, double* vertices, uint8_t* codes,
                    double* cut_vertices, uint8_t* cut_tick_type,
                    double* cut_coords);

/*
 * Final part of a three-part function. Call exactly one time after
 * compute_grid_lines has successfully completed.
 */
int clean_grid_lines_struct(void* struct_ptr);


/*
 * This is the first part of another three-part function.
 * Computes a bounding polygon.
 */
int compute_bounding_polygon(const char* proj_str, double xmin, double xmax,
                             double ymin, double ymax, double atol,
                             double bisection_offset,
                             double minimum_node_distance, void** struct_ptr,
                             size_t* Nvert);

/*
 * Second part of the bounding polygon computation.
 * Saves the bounding polygon.
 */
int save_bounding_polygon(const void* struct_ptr, double* vertices,
                          double* angles);

/*
 * Last part of the bounding polygon computation.
 * Removes the bounding polygon structure.
 */
int clean_bounding_polygon_struct(void* struct_ptr);


}

#endif
