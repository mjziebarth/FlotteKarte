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

#include <../include/griddedinverter.hpp>

using flottekarte::geo_t;
using flottekarte::xy_t;
using flottekarte::GriddedInverter;
using flottekarte::ProjWrapper;

GriddedInverter::GriddedInverter(const ProjWrapper& proj, size_t nlon,
                                 size_t nlat)
{
	const double dlon = 2.0*PI/(nlon+1.0);
	const double dlat = PI/(nlat+1.0);
	for (size_t i=0; i<nlon; ++i){
		for (size_t j=0; j<nlat; ++j){
			const geo_t geo({ (i+0.5) * dlon,
			                  (j+0.5)*dlat - PI/2});
			const xy_t xy(proj.project(geo.lambda, geo.phi));
			tree.insert(std::make_pair(rtpoint(xy.x, xy.y), geo));
		}
	}
}


geo_t GriddedInverter::operator()(const xy_t& xy) const
{
	/* Nearest neighbor: */
	std::vector<rtvalue> result;
	tree.query(boost::geometry::index::nearest(rtpoint(xy.x, xy.y), 1),
	           std::back_inserter(result));
	return result.back().second;
}
