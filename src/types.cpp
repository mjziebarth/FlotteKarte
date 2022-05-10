
#include <../include/types.hpp>
#include <cmath>


double projplot::deg2rad(double deg){
	return projplot::PI / 180.0 * deg;
}

double projplot::rad2deg(double rad){
	return 180.0 / projplot::PI * rad;
}

double projplot::modulo(double a, double b){
	/* True modulo operation (similar to Python's (a % b)).
	 * Implemented here only for positive b (which is what we use).
	 */
	double y = std::fmod(a,b);
	if (y < 0.0)
		return y+b;
	return y;
}
