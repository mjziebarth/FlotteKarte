/*
 * 
 *
 * */

#ifndef PROJPLOT_PYAPI_HPP
#define PROJPLOT_PYAPI_HPP

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

}

#endif
