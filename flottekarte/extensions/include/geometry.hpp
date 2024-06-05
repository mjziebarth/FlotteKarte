/*
 * Geometry.
 *
 * Authors: Malte J. Ziebarth (malte.ziebarth@tum.de)
 *
 * Copyright (C) 2024 Technische Universität München
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

#ifndef FLOTTEKARTE_GEOMETRY_HPP
#define FLOTTEKARTE_GEOMETRY_HPP

#include <../include/types.hpp>

/*
 * Make xy_t available as a boost point:
 */
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/core/cs.hpp>

namespace boost { namespace geometry { namespace traits
{

template<>
struct tag<flottekarte::xy_t> { using type = point_tag; };
template<>
struct dimension<flottekarte::xy_t> : boost::mpl::int_<2> {};
template<>
struct coordinate_type<flottekarte::xy_t> { using type = double; };
template<>
struct coordinate_system<flottekarte::xy_t> {
    using type = boost::geometry::cs::cartesian;
};

template<std::size_t Index>
struct access<flottekarte::xy_t, Index> {
    static_assert(Index < 2, "Out of range");
    using Point = flottekarte::xy_t;
    using CoordinateType = typename coordinate_type<Point>::type;
    static inline CoordinateType get(Point const& p)
    {
        if constexpr (Index == 0)
            return p.x;
        else
            return p.y;
    }
    static inline void set(Point& p, CoordinateType const& value)
    {
        if constexpr (Index == 0)
            p.x = value;
        else
            p.y = value;
    }
};

}}} // namespace boost::geometry::traits

#include <boost/geometry/strategies/strategies.hpp>
#include <boost/geometry/geometries/segment.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/geometry/algorithms/closest_points.hpp>

namespace flottekarte {

/*
 * Line segment
 * ------------
 */
using segment_t = boost::geometry::model::segment<flottekarte::xy_t>;


template<typename index_t>
using SegmentTree
   = boost::geometry::index::rtree<
        segment_t,
        boost::geometry::index::quadratic<16>
>;


template<typename Tree, typename Geometry>
bool within_range(const Tree& tree, const Geometry& geom, double r)
{
    namespace bgi = boost::geometry::index;

    /* Check whether we can query anything: */
    if (tree.size() == 0)
        return false;

    /* Query the nearest segment: */
    typedef typename Tree::value_type value_t;
    std::vector<value_t> nn;
    bgi::query(tree, bgi::nearest(geom, 1), std::back_inserter(nn));

    /* Compute the closest points between the step and the
        * nearest segment: */
    segment_t closest;
    boost::geometry::closest_points(nn[0], geom, closest);

    /* Check whether the distance between the closest points is smaller
        * than r: */
    double dx = closest.first.x - closest.second.x;
    double dy = closest.first.y - closest.second.y;
    return dx*dx + dy*dy <= r*r;
}

}

#endif