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

}

#endif
