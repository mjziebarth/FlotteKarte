/*
 * Gradient computation.
 */

#include <../include/types.hpp>

#ifndef FLOTTEKARTE_GRADIENT_HPP
#define FLOTTEKARTE_GRADIENT_HPP

namespace flottekarte {

enum stencil {
	FORWARD_5POINT
};



template<stencil s>
class Gradient {
public:
	Gradient(const ProjWrapper& proj, const geo_t& coord, double delta=1e-5);

	double gx_east() const;
	double gx_north() const;
	double gy_east() const;
	double gy_north() const;

private:
	double _gx_east;
	double _gx_north;
	double _gy_east;
	double _gy_north;

	void compute_gradients(const ProjWrapper&, const geo_t& coord,
	                       double delta);
};



/*******************************************************************************
 *                          Template implementations.
 ******************************************************************************/

template<>
void Gradient<FORWARD_5POINT>::compute_gradients(const ProjWrapper&,
                                                 const geo_t& coord,
                                                 double delta);


template<stencil s>
double Gradient<s>::gx_east() const
{
	return _gx_east;
}

template<stencil s>
double Gradient<s>::gx_north() const
{
	return _gx_north;
}

template<stencil s>
double Gradient<s>::gy_east() const
{
	return _gy_east;
}

template<stencil s>
double Gradient<s>::gy_north() const
{
	return _gy_north;
}

template<stencil s>
Gradient<s>::Gradient(const ProjWrapper& proj, const geo_t& coord,
                      double delta)
    : _gx_east(0.0), _gx_north(0.0), _gy_east(0.0), _gy_north(0.0)
{
	Gradient<s>::compute_gradients(proj, coord, delta);
}


}

#endif