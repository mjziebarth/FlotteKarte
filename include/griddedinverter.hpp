/*
 * Invert based on gridded forward evaluation. Used mostly to obtain
 * initial conditions for further inversion.
 */

#include <../include/types.hpp>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>

#ifndef FLOTTEKARTE_GRIDDEDINVERTER_HPP
#define FLOTTEKARTE_GRIDDEDINVERTER_HPP

namespace flottekarte {

typedef boost::geometry::cs::cartesian rtcartesian;
typedef boost::geometry::model::point<float, 2, rtcartesian> rtpoint;
typedef std::pair<rtpoint, geo_t> rtvalue;
typedef boost::geometry::index::rstar<16> rtbuild;
typedef boost::geometry::index::rtree<rtvalue, rtbuild> rtree;

class GriddedInverter {
public:
	GriddedInverter(const ProjWrapper& proj, size_t nlon, size_t nlat);

	geo_t operator()(const xy_t& xy) const;

private:
	rtree tree;
};

}

#endif
