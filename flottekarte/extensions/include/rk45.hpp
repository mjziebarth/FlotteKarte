/*
 * Runge-Kutta-Fehlberg integrator rk45
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
 */

#ifndef FLOTTEKARTE_RK45_HPP
#define FLOTTEKARTE_RK45_HPP

#include <array>
#include <optional>
#include <algorithm>
#include <../include/types.hpp>
#include <../include/geometry.hpp>

namespace flottekarte {

namespace rk45 {
/*
 * Coefficients from Table III: "Coefficients for RK4(5), Formula 2"
 * of
 */
constexpr double alpha(uint_fast8_t k)
{
    constexpr std::array<double,6> _alpha = {
        0.0, 1.0 / 4, 3.0 / 8, 12.0 / 13, 1.0, 1.0 / 2
    };
    return _alpha[k];
}

/*
 * beta:
 */
consteval std::array<double,6> beta_l0()
{
    std::array<double, 6> bl0 = {
            0.0,
            1.0 / 4,
            3.0 / 32,
            1932.0 / 2197,
            439.0 / 216,
            -8.0 / 27
    };
    return bl0;
}

consteval std::array<double,6> beta_l1()
{
    std::array<double, 6> bl1 = {
            0.0,
            0.0,
            9.0 / 32,
            -7200.0 / 2197,
            -8.0,
            2.0
    };
    return bl1;
}

consteval std::array<double,6> beta_l2()
{
    std::array<double, 6> bl2 = {
            0.0,
            0.0,
            0.0,
            7296.0 / 2197,
            3680.0 / 513,
            -3544.0 / 2565
    };
    return bl2;
}

consteval std::array<double,6> beta_l3()
{
    std::array<double, 6> bl3 = {
            0.0,
            0.0,
            0.0,
            0.0,
            -845.0 / 4104,
            1859.0 / 4104
    };
    return bl3;
}

consteval std::array<double,6> beta_l4()
{
    std::array<double, 6> bl4 = {
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            -11.0 / 40
    };
    return bl4;
}


constexpr double beta(uint_fast8_t k, uint_fast8_t l)
{
    constexpr std::array<std::array<double,6>, 5> _beta = {
        beta_l0(),
        beta_l1(),
        beta_l2(),
        beta_l3(),
        beta_l4()
    };
    return _beta.at(l).at(k);
}

constexpr double c(uint_fast8_t k)
{
    if (k == 0)
        return 25.0 / 216;
    else if (k == 1)
        return 0.0;
    else if (k == 2)
        return 1408.0 / 2565;
    else if (k == 3)
        return 2197.0 / 4104;
    else if (k == 4)
        return - 1.0 / 5;
    return 0.0;
}

constexpr double c_hat(uint_fast8_t k)
{
    if (k == 0)
        return 16.0 / 135;
    else if (k == 1)
        return 0.0;
    else if (k == 2)
        return 6656.0 / 12825;
    else if (k == 3)
        return 28561.0 / 56430;
    else if (k == 4)
        return - 9.0 / 50;
    else if (k == 5)
        return 2.0 / 55;
    return 0.0;
}

constexpr double d_te(uint_fast8_t k)
{
    if (k == 0)
        return -1.0 / 360;
    else if (k == 1)
        return 128.0 / 4275;
    else if (k == 2)
        return 2187.0 / 75240;
    else if (k == 3)
        return -1.0 / 50;
    else if (k == 4)
        return -2.0 / 55;

    return 0.0;
}

} // namespace rk45

template<
    int direction,
    typename fun_t,
    typename exit_cond_t,
    typename step_index_t,
    bool prevent_self_intersections=true
>
void
forward_rk45(
    fun_t function,
    exit_cond_t exit,
    double x0,
    double y0,
    double xmin,
    double xmax,
    double ymin,
    double ymax,
    step_index_t Nmax,
    double epsilon,
    double r,
    double ds_min,
    SegmentTree<step_index_t>& stree,
    path_xy_t& history
)
{
    namespace bgi = boost::geometry::index;

    auto in_bounds = [xmin, xmax, ymin, ymax](const xy_t& xy) -> bool
    {
        return xy.x >= xmin && xy.x <= xmax && xy.y >= ymin && xy.y <= ymax;
    };
    auto norm = [](const xy_t& xy) -> double
    {
        return std::sqrt(xy.x*xy.x + xy.y*xy.y);
    };

    auto r_self_intersection = [r, &stree](const segment_t& s) -> bool
    {
        if constexpr (!prevent_self_intersections)
            return false;

//        std::cerr << "Actually implement self-intersection test!\n";
        return false;
        throw std::runtime_error("Actually implement self-intersection test!");
        return within_range(stree, s, r);
    };

    xy_t xy_i(x0, y0);

    /* Ensure we start inside bounds: */
    if (!in_bounds(xy_i))
        return;

    history.push_back(xy_i);

    constexpr double H_INITIAL = 1e-3;
    double h_abs = H_INITIAL;
    std::array<xy_t,6> k;

    /* Here we sum the distance to the last point: */
    double ds = 0.0;

    for (size_t i=0; i<Nmax; ++i)
    {
        /* Directional step (forward or backward): */
        double h = direction * h_abs;

        /* Compute the direction vector: */
        xy_t v(function(xy_i));

        /* Compute the k: */
        k[0] = h * function(xy_i);
        for (uint_fast8_t j=1; j<6; ++j){
            xy_t xy_j(xy_i);
            for (uint_fast8_t l=0; l < j; ++l){
                xy_j += rk45::beta(j, l) * k[l];
            }
            k[j] = h * function(xy_j);
        }

        /* Estimate of the truncation error: */
        xy_t x_te(rk45::d_te(0) * k[0]);
        for (uint_fast8_t j=1; j<6; ++j)
            x_te += rk45::d_te(j) * k[j];
        double TE = norm(x_te);

        /* New step width: */
        h_abs = 0.9 * h_abs * std::pow(epsilon / TE, 1.0/5);
        if (std::isinf(h_abs) || std::isnan(h_abs))
            h_abs = H_INITIAL;

        /* If the approximation error is too large, need to redo this
         * step with the shorter step size: */
        if (TE > epsilon)
            continue;

        /* Accepted the new point: */
        xy_t xy_next = xy_i + rk45::c(0) * k[0];
        for (uint_fast8_t j=1; j<6; ++j)
            xy_next += rk45::c(j) * k[j];

        /* Check whether the new point intersects any existing segment: */
        segment_t seg_i(xy_i, xy_next);
        if (r_self_intersection(seg_i)){
            /* We went out of bounds. Crop this segment to the bounding
             * box: */

            /* Stop the integration. */
            break;
        }

        if (exit(xy_next))
            /* Triggered the exit condition */
            break;


        if (in_bounds(xy_next)){
            /* Add to the distance traveled: */
            ds += norm(xy_next - xy_i);

            /* Add segment to history if distance is sufficient: */
            if (ds >= ds_min){
                ds = 0.0;
                history.push_back(xy_next);
            }

            if constexpr (prevent_self_intersections){
                /* Add to segment tree: */
                stree.insert(
                    segment_t(xy_i, xy_next)
                );
            }

            /* Update location: */
            xy_i = xy_next;
        } else {
            /* We went out of bounds. Crop this segment to the bounding
             * box: */

            /* Stop the integration: */
            break;
        }

    }

    return;
}

template<typename fun_t, typename exit_cond_t, typename step_index_t>
path_xy_t
forward_backward_rk45(
    fun_t function,
    exit_cond_t exit,
    double x0,
    double y0,
    double xmin,
    double xmax,
    double ymin,
    double ymax,
    step_index_t Nmax,
    double epsilon,
    double r,
    double ds_min
)
{

    SegmentTree<step_index_t> stree;

    path_xy_t path;
    forward_rk45<-1>(
        function, exit, x0, y0, xmin, xmax, ymin, ymax,
        Nmax, epsilon, r, ds_min, stree, path
    );

    if (path.size() > 0){
        /* Start at the end point: */
        std::reverse(path.begin(), path.end());

        /* Do not keep the start point (it will be added again now): */
        path.resize(path.size() - 1);
    }
    /* Append the forward integration: */
    forward_rk45<1>(
        function, exit, x0, y0, xmin, xmax, ymin, ymax,
        Nmax, epsilon, r, ds_min, stree, path
    );

    return path;
}

}


#endif