/*
 * Gradient computation.
 *
 * Authors: Malte J. Ziebarth (ziebarth@gfz-potsdam.de)
 *
 * Copyright (C) 2022 Deutsches GeoForschungsZentrum Potsdam
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