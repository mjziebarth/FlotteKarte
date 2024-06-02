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

namespace flottekarte {

struct segment_t;
struct circle_t;

/*
 * Axis-aligned bounding box
 * -------------------------
 */
struct bbox_t {
    double xmin;
    double xmax;
    double ymin;
    double ymax;

    bbox_t(double xmin, double xmax, double ymin, double ymax);

    bbox_t(const xy_t& p0, const xy_t& p1);

    bbox_t(const bbox_t& b0, const bbox_t& b1);

    bbox_t(const segment_t& s);

    bbox_t(const circle_t& c);

    double area() const;
};


/*
 * Line segment
 * ------------
 */
struct segment_t {
    xy_t p0;
    xy_t p1;

    segment_t() = default;
    segment_t(const xy_t& p0, const xy_t& p1);

    bbox_t bbox() const;

    bool operator==(const segment_t& other) const;
};


/*
 * Circle
 * ------
 */
struct circle_t {
    xy_t center;
    double r = 0.0;

    circle_t() = default;
    circle_t(const xy_t& center, double r);

    bbox_t bbox() const;

    bool operator==(const circle_t& other) const;
};


/*
 * Distances:
 */
double distance2(const xy_t& p, const segment_t& s);
double distance(const xy_t& p, const segment_t& s);


/*
 * Spatial relations:
 */
bool contains(const bbox_t& b0, const bbox_t& b1);

bool intersects(const bbox_t& b0, const bbox_t& b1);
bool intersects(const segment_t& s0, const segment_t& s1);

bool intersects(const circle_t& c, const segment_t& s);
bool intersects(const segment_t& s, const circle_t& c);


}

#endif