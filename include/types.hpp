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

struct geo_degrees_t {
	double lon;
	double lat;
};

enum axis_t {
	AX_BOT=0, AX_TOP=1, AX_LEFT=2, AX_RIGHT=3
};

enum tick_t {
	TICK_LON=0, TICK_LAT=1
};


double deg2rad(double deg);
double rad2deg(double rad);

double modulo(double a, double b);

}

#endif