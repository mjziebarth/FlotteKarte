/*
 *
 */

#include <../include/griddedinverter.hpp>

using projplot::geo_t;
using projplot::xy_t;
using projplot::GriddedInverter;
using projplot::ProjWrapper;

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
