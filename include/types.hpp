/*
 * Some common types.
 */

#ifndef PROJPLOT_TYPES_HPP
#define PROJPLOT_TYPES_HPP

namespace projplot {

constexpr double PI = 3.14159265358979323846;

struct xy_t {
	double x;
	double y;
};

struct geo_t {
	double lambda;
	double phi;
};

double deg2rad(double deg);
double rad2deg(double rad);

double modulo(double a, double b);

}

#endif
