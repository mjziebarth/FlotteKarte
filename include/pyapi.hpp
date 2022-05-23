/*
 * 
 *
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

int compute_axes_ticks(const char* proj_str, double xmin, double xmax,
                       double ymin, double ymax, double tick_spacing_degree,
                       unsigned int max_ticks_per_axis,
                       unsigned char which_bot, unsigned char which_top,
                       unsigned char which_left, unsigned char which_right,
                       double* bot_ticks, double* top_ticks,
                       double* left_ticks, double* right_ticks,
                       unsigned int* Nticks);

/*
 * This is the first of a three-part function. Computes grid lines according to
 * settings (projection, x- & ylim, tick spacing, bisection offset,
 * minimum distance between path-adjacent nodes) and returns 1) the length
 * of the resulting path (Npath) and 2) a pointer to the structure holding
 * all the required information (struct_ptr).
 * It should to be followed, after allocating two suitable numpy arrays,
 * by a call to save_grid_lines, which transfers the path contained in
 * struct_ptr and frees the allocated space.
 * It then has to be followed by a call to clean_grid_lines_struct. If this
 * last call is omitted, a memory leak will occur. If this is called a second
 * time, deallocated memory will be accessed.
 */
int compute_grid_lines(const char* proj_str, double xmin, double xmax,
                       double ymin, double ymax, int tick_spacing_degree,
                       double bisection_offset, double minimum_node_distance,
                       double max_lat, double cut_at_angle_degrees,
                       void** struct_ptr, size_t* Npath);

/*
 * Second part of a three-part function.
 */
int save_grid_lines(const void* struct_ptr, double* vertices, uint8_t* codes);

/*
 * Final part of a three-part function. Call exactly one time after
 * compute_grid_lines has successfully completed.
 */
int clean_grid_lines_struct(void* struct_ptr);

}

#endif
