
#include <../include/types.hpp>
#include <cmath>

using flottekarte::xy_t;
using flottekarte::geo_t;
using flottekarte::geo_degrees_t;
using flottekarte::deg2rad;



double xy_t::distance(const xy_t& other) const
{
	return std::sqrt(distance_squared(other));
}

double xy_t::distance_squared(const xy_t& other) const
{
	const double dx = x - other.x;
	const double dy = y - other.y;
	return dx*dx + dy*dy;
}


bool xy_t::operator==(const xy_t& other) const
{
	return (x == other.x) && (y == other.y);
}


geo_degrees_t::geo_degrees_t(double lon, double lat) : lon(lon), lat(lat)
{
}


geo_t geo_degrees_t::to_radians() const
{
	return {deg2rad(lon), deg2rad(lat)};
}


double flottekarte::deg2rad(double deg){
	return flottekarte::PI / 180.0 * deg;
}

double flottekarte::rad2deg(double rad){
	return 180.0 / flottekarte::PI * rad;
}

double flottekarte::modulo(double a, double b){
	/* True modulo operation (similar to Python's (a % b)).
	 * Implemented here only for positive b (which is what we use).
	 */
	double y = std::fmod(a,b);
	if (y < 0.0)
		return y+b;
	return y;
}
