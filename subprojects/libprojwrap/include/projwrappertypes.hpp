/*
 *
 */

#ifndef PROJWRAPPER_PROJWRAPPERTYPES_HPP
#define PROJWRAPPER_PROJWRAPPERTYPES_HPP

namespace projwrapper {

constexpr double PI = 3.14159265358979323846;

struct xy_t {
	double x;
	double y;

	double distance(const xy_t& other) const;

	double distance_squared(const xy_t& other) const;

	bool operator==(const xy_t& other) const;
};

struct geo_t {
	double lambda;
	double phi;
};

struct geo_degrees_t {
	double lon;
	double lat;

	geo_degrees_t(double lon, double lat);

	geo_t to_radians() const;
};


double deg2rad(double deg);
double rad2deg(double rad);

double modulo(double a, double b);

}

#endif