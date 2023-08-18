/*
 * Invert based on gridded forward evaluation. Used mostly to obtain
 * initial conditions for further inversion.
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
