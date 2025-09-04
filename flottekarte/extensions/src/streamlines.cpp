/*
 * Generation of streamline polygons.
 *
 * Authors: Malte J. Ziebarth (malte.ziebarth@tum.de)
 *
 * Copyright (C) 2024-2025 Technische Universität München
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


static path_xy_t line_buffer(
    const path_xy_t& line,
    const std::vector<double>& widths
)
{
    /* Need at least one point in the middle: */
    if (line.size() < 3)
        return path_xy_t();

    if (widths.size()+2 != line.size())
        throw std::runtime_error("Internal vector corruption in line_buffer "
            "arguments."
        );

    /* In two vectors, we place the left and the right offset paths.
     * Once we have traversed the whole 'line', we merge them. */
    path_xy_t left(0), right(0);
    left.reserve(line.size() - 2);
    right.reserve(2*line.size() - 2);

    /* Start the right path with the first point of the path: */
    auto it = line.cbegin();
    xy_t x_i(*it);
    right.emplace_back(x_i);
    ++it;
    xy_t last_direction(*it - x_i);

    /* Now iterate through all interior points: */
    const auto end = line.cend() - 1;
    auto wit = widths.cbegin();

    while (it != end)
    {
        /* Get the target width: */
        double width =*wit;
        if (std::isinf(width))
            throw std::runtime_error("width is inf!");

        /* Get the next direction by advancing through the line string: */
        xy_t x_i(*it);
        ++it;
        ++wit;
        xy_t next_direction(*it - x_i);

        /* Vector that points left from the current average direction.
         * Norm it to the correct distance. */
        xy_t avg_dir(last_direction + next_direction);
        xy_t vl(-avg_dir.y, avg_dir.x);
        vl *= width / vl.norm();

        /* Add the two points: */
        left.push_back(x_i + vl);
        right.push_back(x_i - vl);

        /* Advance: */
        last_direction = next_direction;
    }

    /* We started with the first point in right, add the final point as well: */
    right.push_back(*end);

    /* Append the way back on the left side: */
    std::reverse(left.begin(), left.end());
    right.insert(right.end(), left.begin(), left.end());

    return right;
}

/*
 *
 * Streamlines main method:
 *
 */
std::vector<path_xy_t>
streamlines(
    double xmin,
    double xmax,
    size_t nx,
    double ymin,
    double ymax,
    size_t ny,
    const double* z,
    double r,
    double ds_min,
    double width_scale,
    double epsilon,
    width_mode_t width_mode
)
{
    if (width_scale <= 0.0)
        throw std::runtime_error("width_scale has to be positive!");

    constexpr uint16_t STARTING_POINTS = 512;
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


        /* The vector field. */
        auto direction
        = [&interp](const xy_t& p) -> xy_t
        {
            double alpha = interp(p)[0];
            return xy_t(std::sin(alpha), std::cos(alpha));
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

    /* Now we need to transform the paths into polygons of appropriate width
     * around the integration path.
     * Start with a function that computes the width as a function of location:
     */
    auto compute_width
    = [width_mode, r, width_scale, &interp]
    (const xy_t& p) -> double
    {
        std::array<double,3> vf_p = interp(p);
        double p0 = vf_p[1];
        double p1 = vf_p[2];
        double raw_result = 0.0;
        if (width_mode == DIFFERENCE)
            raw_result = std::abs(p0 - p1);
        else if (width_mode == SUM)
            raw_result = std::abs(p0 + p1);
        else if (width_mode == CONSTANT)
            raw_result = 1.0;
        else
            throw std::runtime_error("Unknown width_mode.");

        // Here we hard-code a 10% margin to filling the whole plane:
        return 0.95 * r * raw_result / width_scale;
    };
    std::vector<double> widths;
    std::vector<size_t> nans;
    std::vector<path_xy_t> split_off;
    std::vector<size_t> to_delete;
    for (size_t p=0; p<result.size(); ++p){
        path_xy_t& path = result[p];
        /* Compute the corresponding widths: */
        widths.resize(path.size()-2);
        for (size_t i=1; i<path.size()-1; ++i)
            widths[i-1] = compute_width(path[i]);

        /* Check if any of them are nan or inf: */
        nans.clear();
        for (size_t i=0; i<widths.size(); ++i){
            if (std::isnan(widths[i]) || std::isinf(widths[i]))
                nans.push_back(i+1);
        }

        /* If any nan, we have to split. Otherwise, we can just invoke the
         * buffer once: */
        if (nans.size() == 0)
            path = line_buffer(path, widths);
        else {
            /* We need to split the original path at all NaN occurrences. */
            size_t i0 = 0;
            path_xy_t path_original;
            path_original.swap(path);
            size_t n=0;
            /* Add, artificially, the end index to nan so that we don't have
             * to treat the last index interval separately: */
            nans.push_back(path_original.size());
            for (size_t i1 : nans){
                if (i1 - i0 >= Nmin){
                    std::vector<double> wtmp(widths.cbegin()+i0, widths.cbegin()+i1-2);
                    path_xy_t ptmp(
                        path_original.cbegin()+i0,
                        path_original.cbegin()+i1
                    );
                    if (n == 0)
                        path = line_buffer(ptmp, wtmp);
                    else
                        split_off.emplace_back(line_buffer(ptmp, wtmp));

                    ++n;
                }
                i0 = i1;
            }
            /* Everything might be nan. Then delete this path: */
            if (n == 0)
                to_delete.push_back(p);
        }
    }

    /* Add the split off paths: */
    if (split_off.size() > 0){
        size_t n = result.size();
        result.resize(n + split_off.size());
        auto to = result.begin() + n;
        for (auto from = split_off.begin(); from != split_off.end(); ++from)
        {
            from->swap(*to);
            ++to;
        }
    }


    /* Delete paths.
     * Create new vector, swap all non-delete paths, and swap the
     * new vector back to results. */
    if (to_delete.size() > 0){
        /*
         * Did we delete all streamlines?
         */
        if (to_delete.size() == result.size()){
            result.clear();
            return result;
        }
        std::vector<path_xy_t> tmp(result.size() - to_delete.size());
        auto from = result.begin();
        auto to = tmp.begin();
        auto dit = to_delete.cbegin();
        size_t i=0;
        for (; i < result.size(); ++i){
            if (i == *dit){
                /* We are currently at a node we do not want to delete. */
                ++dit;
            } else {
                if (to == tmp.cend())
                    throw std::runtime_error("to == tmp.cend() (loc 1)");
                from->swap(*to);
                ++to;
            }
            ++from;
            if (dit == to_delete.cend())
                break;
        }

        /* After the last one: */
        for (; from != result.cend(); ++from)
        {
            if (to == tmp.cend())
                throw std::runtime_error("to == tmp.cend()");
            from->swap(*to);
            ++to;
        }

        result.swap(tmp);
    }

    return result;
}


}