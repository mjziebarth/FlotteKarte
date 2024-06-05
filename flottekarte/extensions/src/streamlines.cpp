/*
 * Generation of streamline polygons.
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
 * */

#include <../include/interpolate.hpp>
#include <../include/streamlines.hpp>
#include <../include/rk45.hpp>

#include <boost/random/sobol.hpp>
#include <boost/random/uniform_real_distribution.hpp>

#include <algorithm>

namespace flottekarte {

std::vector<path_xy_t>
streamlines(
    double xmin, double xmax, size_t nx,
    double ymin, double ymax, size_t ny,
    const double* z,
    double r, double ds_min,
    double epsilon, bool azimuth_is_full_circle
)
{

    constexpr uint16_t STARTING_POINTS = 100;
    constexpr uint32_t Nmin = 5;

    const uint32_t Nmax = 100000;

    /* Set the precision relative to the spatial extent: */
    epsilon *= std::max(xmax - xmin, ymax-ymin);

    /* Setup the interplator: */
    LinearGrid2DInterpolator<false,3> interp(xmin, xmax, nx, ymin, ymax, ny, z);

    /* The resulting path: */
    std::vector<path_xy_t> result;

    /* Segment tree to query for intersection with previous segments: */
    SegmentTree<uint32_t> segment_tree;

    /* Get the starting points using a Sobol sequence: */
    std::vector<xy_t> start_points(STARTING_POINTS);
    boost::random::sobol qrng(2);
    boost::random::uniform_real_distribution<double> distx(xmin, xmax);
    boost::random::uniform_real_distribution<double> disty(ymin, ymax);
    for (size_t i=0; i<STARTING_POINTS; ++i)
    {
        start_points[i].x = distx(qrng);
        start_points[i].y = disty(qrng);
    }


    auto early_exit = [r, &segment_tree](const xy_t& p) -> bool
    {
        return within_range(segment_tree, p, r);
    };


    for (uint16_t i=0; i<STARTING_POINTS; ++i)
    {
        /* Set the starting point: */
        xy_t start = start_points[i];


        /* The vector field.
         * Here we might have to perform some tricks.
         * If the azimuth is only half range (180°), we might have 180°
         * jump discontinuities in the orientation. */
        std::array<xy_t, 12> direction_memory;
        uint_fast8_t mem_id = 0;
        xy_t last_direction;
        {
            double last_alpha = interp(start)[0];
            last_direction.x = std::sin(last_alpha);
            last_direction.y = std::cos(last_alpha);
        }
        auto direction
        = [azimuth_is_full_circle, &interp, &direction_memory, &mem_id]
        (const xy_t& p) -> xy_t
        {
            std::array<double,3> vf_p = interp(p);
            double alpha = vf_p[0];
            xy_t dir(std::sin(alpha), std::cos(alpha));
            /* If the new direction is rather opposite the old direction,
             * and the azimuth is not given in full range (i.e. we might have
             * jump discontinuities), flip the direction here: */
            if (!azimuth_is_full_circle){
                xy_t avg_dir(std::accumulate(
                    direction_memory.cbegin(), direction_memory.cend(),
                    xy_t(0.0, 0.0)
                ));
                if (avg_dir.dot(dir) < 0.0){
                    dir *= -1.0;
                }
            }
            direction_memory[mem_id++] = dir;
            if (mem_id == 12)
                mem_id = 0;
            return dir;
        };

        /* Integrate: */
        path_xy_t path = forward_backward_rk45(
            direction,
            early_exit,
            start.x, start.y,
            xmin, xmax, ymin, ymax, Nmax,
            epsilon, r, ds_min
        );

        /* Take only paths of minimum size: */
        if (path.size() < Nmin)
            continue;

        /* Add the segments to the tree: */
        auto it = path.cbegin();
        xy_t last = *it;
        for (++it; it != path.cend(); ++it)
        {
            segment_tree.insert(segment_t(last, *it));
            last = *it;
        }

        /* Append the path to the result: */
        result.emplace_back();
        result.back().swap(path);
    }

    return result;
}


}