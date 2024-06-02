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

#include <iostream>
#include <cmath>
#include <../include/geometry.hpp>

namespace flottekarte {


bbox_t::bbox_t(double xmin, double xmax, double ymin, double ymax)
   : xmin(xmin), xmax(xmax), ymin(ymin), ymax(ymax)
{}


bbox_t::bbox_t(const xy_t& p0, const xy_t& p1)
{
    if (p0.x <= p1.x){
        xmin = p0.x;
        xmax = p1.x;
    } else {
        xmin = p1.x;
        xmax = p0.x;
    }
    if (p0.y <= p1.y){
        ymin = p0.y;
        ymax = p1.y;
    } else {
        ymin = p1.y;
        ymax = p0.y;
    }
}

bbox_t::bbox_t(const segment_t& s)
{
    if (s.p0.x <= s.p1.x){
        xmin = s.p0.x;
        xmax = s.p1.x;
    } else {
        xmin = s.p1.x;
        xmax = s.p0.x;
    }
    if (s.p0.y <= s.p1.y){
        ymin = s.p0.y;
        ymax = s.p1.y;
    } else {
        ymin = s.p1.y;
        ymax = s.p0.y;
    }
}

bbox_t::bbox_t(const bbox_t& b0, const bbox_t& b1)
{
    if (b0.xmin <= b1.xmin)
        xmin = b0.xmin;
    else
        xmin = b1.xmin;
    if (b0.xmax >= b1.xmax)
        xmax = b0.xmax;
    else
        xmax = b1.xmax;
    if (b0.ymin <= b1.ymin)
        ymin = b0.ymin;
    else
        ymin = b1.ymin;
    if (b0.ymax >= b1.ymax)
        ymax = b0.ymax;
    else
        ymax = b1.ymax;

}

bbox_t::bbox_t(const circle_t& c)
   : xmin(c.center.x - c.r), xmax(c.center.x + c.r), ymin(c.center.y - c.r),
     ymax(c.center.y + c.r)
{}


double bbox_t::area() const
{
    return (xmax - xmin) * (ymax - ymin);
}


/*
 * Segment
 */
segment_t::segment_t(const xy_t& p0, const xy_t& p1) : p0(p0), p1(p1)
{
}


bool segment_t::operator==(const segment_t& other) const
{
    return p0 == other.p0 && p1 == other.p1;
}


bbox_t segment_t::bbox() const
{
    return bbox_t(p0, p1);
}

/*
 * Circle
 */

circle_t::circle_t(const xy_t& center, double r) : center(center), r(r)
{}

bbox_t circle_t::bbox() const
{
    return bbox_t(center.x-r, center.x+r, center.y-r, center.y+r);
}

bool circle_t::operator==(const circle_t& c) const
{
    return center == c.center && r == c.r;
}



/*
 * Distances:
 */
double distance2(const xy_t& p, const segment_t& s)
{
    /* This follows
     *    Morrison, Jack C. (1995): "Distance from a point to a line". Graphics
     *    Gems II.
     * with:
     *    P = p
     *    A = s.p0
     *    B = s.p1
     */
    long double Ax = s.p0.x;
    long double Ay = s.p0.y;
    long double Bx = s.p1.x;
    long double By = s.p1.y;
    long double PAx = p.x - Ax;
    long double PAy = p.y - Ay;
    long double BAx = Bx - Ax;
    long double BAy = By - Ay;
    /* Check if p beyond s.p0: */
    long double t = PAx * BAx + PAy * BAy;
    if (t < 0)
        return PAx * PAx + PAy * PAy;
    long double BPx = Bx - p.x;
    long double BPy = By - p.y;
    /* Check if p beyond s.p1: */
    t = BPx * BAx + BPy * BAy;
    if (t < 0)
        return BPx * BPx + BPy * BPy;

    /* Compute the distance for p in between: */
    double a2 = PAy * BAx - PAx * BAy;
    double BAx2 = BAx*BAx;
    double BAy2 = BAy*BAy;
    double d1_2 = a2 * a2 / (BAx2 + BAy2);

    return d1_2;
}

double distance(const xy_t& p, const segment_t& s)
{
    return std::sqrt(distance2(p,s));
}


/*
 * Spatial set relations:
 */


bool contains(const bbox_t& b0, const bbox_t& b1)
{
    return b1.xmin >= b0.xmin && b1.xmax <= b0.xmax
        && b1.ymin >= b0.ymin && b1.ymax <= b0.ymax;
}


bool intersects(const bbox_t& b0, const bbox_t& b1)
{
    bool not_intersect = (
        b0.xmin > b1.xmax || b0.xmax < b1.xmin ||
        b0.ymin > b1.ymax || b0.ymax < b1.ymin
    );
    return !not_intersect;
}

bool intersects(const segment_t& s0, const segment_t& s1)
{
    /* Follows Antonio, Franklin (1992): Faster Line Segment Intersection */
    /* P1 = s0.p0,
     * P2 = s0.p1,
     * P3 = s1.p0,
     * P4 = s1.p1
     * A = P2 - P1 = s0.p1 - s0.p0,
     * B = P3 - P4 = s1.p0 - s1.p1,
     * C = P1 - P3 = s0.p0 - s1.p0 */
    const long double P1x = s0.p0.x;
    const long double P1y = s0.p0.y;
    const long double P2x = s0.p1.x;
    const long double P2y = s0.p1.y;
    const long double P3x = s1.p0.x;
    const long double P3y = s1.p0.y;
    const long double P4x = s1.p1.x;
    const long double P4y = s1.p1.y;
    const long double Ax = P2x - P1x;
    const long double Ay = P2y - P1y;
    const long double Bx = P3x - P4x;
    const long double By = P3y - P4y;
    const long double Cx = P1x - P3x;
    const long double Cy = P1y - P3y;
    /* Degenerate cases first: */
    if (Ax == 0.0 && Ay == 0.0)
        return false;
    if (Bx == 0.0 && By == 0.0)
        return false;
    if ((Ax == 0.0 && Bx == 0.0) || (Ay == 0.0 && By == 0.0))
        return intersects(bbox_t(s0), bbox_t(s1));

    const double denom = (Ay * Bx - Ax * By);
    if (denom > 0){
        const double alpha_num = (By * Cx - Bx * Cy);
        if (alpha_num < 0 || alpha_num > denom)
            return false;

        const double beta_num  = (Ax * Cy - Ay * Cx);
        if (beta_num < 0 || beta_num > denom)
            return false;
    } else {
        const double alpha_num = (By * Cx - Bx * Cy);
        if (alpha_num > 0 || alpha_num < denom)
            return false;

        const double beta_num  = (Ax * Cy - Ay * Cx);
        if (beta_num > 0 || beta_num < denom)
            return false;
    }
    return true;
}

/*
 * Circle-segment intersection
 */
bool intersects(const circle_t& c, const segment_t& s)
{
    /* Circle with center x and intersects with a line segment iff the line
     * segment's closest point to x is at a distance less than or equal to the
     * circle's radius r. */
    double r2 = c.r * c.r;
    return distance2(c.center, s) <= r2;
}

bool intersects(const segment_t& s, const circle_t& c)
{
    return intersects(c, s);
}


} // end namespace