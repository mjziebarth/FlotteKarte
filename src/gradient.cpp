/*
 *
 */

#include <../include/gradient.hpp>
#include <cmath>

using flottekarte::xy_t;
using flottekarte::geo_t;
using flottekarte::Gradient;
using flottekarte::ProjWrapper;
using flottekarte::FORWARD_5POINT;

#include <iostream>

template<>
void Gradient<FORWARD_5POINT>::compute_gradients(const ProjWrapper& proj,
                                            const geo_t& coord, double delta)
{
	/* Use five-point forward stencil to compute 4th order forward
	 * finite difference. */
	const double a = proj.a();
	const double f = proj.f();
	const double e2 = f*(2.0 - f);

	/* Make sure that we do not exceed the latitude bounds: */
	const double sign_lat = (coord.phi + 4 * delta > PI / 2) ? -1.0 : 1.0;
	const double sign_lon = (coord.lambda + 4 * delta > PI) ? -1.0 : 1.0;

	/* The reference: */
	const xy_t xy(proj.project(coord));

	/* Longitude stencil: */
	const double idel = 1.0 / delta;
	const double& la = coord.lambda;
	const double& ph = coord.phi;
	const xy_t xy_dlon1(proj.project(geo_t({la +   sign_lon*delta, ph})));
	const xy_t xy_dlon2(proj.project(geo_t({la + 2*sign_lon*delta, ph})));
	const xy_t xy_dlon3(proj.project(geo_t({la + 3*sign_lon*delta, ph})));
	const xy_t xy_dlon4(proj.project(geo_t({la + 4*sign_lon*delta, ph})));
	_gx_east = (-0.25 * xy_dlon4.x + 4.0/3.0 * xy_dlon3.x - 3.0 * xy_dlon2.x
	            + 4.0 * xy_dlon1.x - 25.0/12.0 * xy.x) * idel;
	_gy_east = (-0.25 * xy_dlon4.y + 4.0/3.0 * xy_dlon3.y - 3.0 * xy_dlon2.y
	            + 4.0 * xy_dlon1.y - 25.0/12.0 * xy.y) * idel;

	/* Latitude stencil: */
	const xy_t xy_dlat1(proj.project(geo_t({la, ph +   sign_lat*delta})));
	const xy_t xy_dlat2(proj.project(geo_t({la, ph + 2*sign_lat*delta})));
	const xy_t xy_dlat3(proj.project(geo_t({la, ph + 3*sign_lat*delta})));
	const xy_t xy_dlat4(proj.project(geo_t({la, ph + 4*sign_lat*delta})));
	_gx_north = (-0.25 * xy_dlat4.x + 4.0/3.0 * xy_dlat3.x - 3.0 * xy_dlat2.x
	             + 4.0 * xy_dlat1.x - 25.0/12.0 * xy.x) * idel;
	_gy_north = (-0.25 * xy_dlat4.y + 4.0/3.0 * xy_dlat3.y - 3.0 * xy_dlat2.y
	             + 4.0 * xy_dlat1.y - 25.0/12.0 * xy.y) * idel;

	/* Now handle the transformation of lon/lat onto the ellipsoid surface.
	 * Also correct the sign for backwards finite difference that we may have
	 * introduced.
	 */
	const double sp = std::sin(coord.phi);
	const double B = std::sqrt(1.0 - e2*sp*sp);
	const double transform_east = B / (a * std::cos(coord.phi)) * sign_lon;
	_gx_east *= transform_east;
	_gy_east *= transform_east;
	const double transform_north = B*B*B / (a * (1.0 - e2)) * sign_lat;
	_gx_north *= transform_north;
	_gy_north *= transform_north;
}
