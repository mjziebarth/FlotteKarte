/*
 * Path processing facilities.
 *
 * Authors: Malte J. Ziebarth (ziebarth@gfz-potsdam.de)
 *
 * Copyright (C) 2022 Deutsches GeoForschungsZentrum Potsdam,
 *                    Malte J. Ziebarth
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

#include <../include/paths.hpp>



#include <optional>
#include <stack>
#include <iostream>
#include <cmath>
#include <set>


using flottekarte::xy_t;
using flottekarte::geo_t;
using flottekarte::geo_degrees_t;
using flottekarte::path_xy_t;
using flottekarte::refined_path_t;
using flottekarte::deg2rad;
using flottekarte::ProjWrapper;
using flottekarte::crop_and_refine;



xy_t flottekarte::boundary_intersection(const xy_t& out, const xy_t& in,
                                        double xmin, double xmax, double ymin,
                                        double ymax)
{
	/* These functions interpolates the coordinates of outside and inside
	 * points once the crossed boundary has been identified and the distances
	 * dz_out and dz_in from that boundary have been computed.
	 */
	auto interpolate_x = [&](double dz_out, double dz_in) -> double {
		double wtot = dz_out + dz_in;
		double wout = dz_in / wtot;
		double win = dz_out / wtot;
		return wout * out.x + win * in.x;
	};
	auto interpolate_y = [&](double dz_out, double dz_in) -> double {
		double wtot = dz_out + dz_in;
		double wout = dz_in / wtot;
		double win = dz_out / wtot;
		return wout * out.y + win * in.y;
	};

	/* Determine the intersection of a line (out,in) with the bounding box
	 * given by (xmin,xmax,ymin,ymax) */
	if (out.y > ymax){
		/* Intersection with the y axis: */
		const double xyint = interpolate_x(out.y - ymax, ymax - in.y);
		if (out.x < xmin){
			/* Intersection with the x axis: */
			double yxint = interpolate_y(xmin - out.x, in.x - xmin);
			if (yxint < ymax){
				/* Intersection with x axis leaves the bounding box: */
				return {xmin, yxint};
			}
			/* Intersection with the y axis leaves the bounding box
			 * (or corner intersection): */
			return {xyint, ymax};

		} else if (out.x > xmax){
			/* Intersection with the x axis: */
			double yxint = interpolate_y(out.x - xmax, xmax - in.x);
			if (yxint < ymax){
				/* Intersection with x axis leaves the bounding box: */
				return {xmax, yxint};
			}
			/* Intersection with the y axis leaves the bounding box
			 * (or corner intersection): */
			return {xyint, ymax};

		}
		/* The simple case that only y > ymax. */
		return {xyint, ymax};

	} else if (out.y < ymin){
		const double xyint = interpolate_x(ymin - out.y, in.y - ymin);
		if (out.x < xmin){
			/* Intersection with the x axis: */
			double yxint = interpolate_y(xmin - out.x, in.x - xmin);
			if (yxint > ymin){
				/* Intersection with x axis leaves the bounding box: */
				return {xmin, yxint};
			}
			/* Intersection with the y axis leaves the bounding box
			 * (or corner intersection): */
			return {xyint, ymin};

		} else if (out.x > xmax){
			/* Intersection with the x axis: */
			double yxint = interpolate_y(out.x - xmax, xmax - in.x);
			if (yxint > ymin){
				/* Intersection with x axis leaves the bounding box: */
				return {xmax, yxint};
			}
			/* Intersection with the y axis leaves the bounding box
			 * (or corner intersection): */
			return {xyint, ymin};

		}
		/* The simple case that only y < ymin. */
		return {xyint, ymin};

	} else if (out.x > xmax){
		/* The simple case that only x > xmax. */
		return {xmax, interpolate_y(out.x - xmax, xmax - in.x)};

	}
	/* The simple case that only x < xmin. */
	return {xmin, interpolate_y(xmin - out.x, in.x - xmin)};

}


static double distance_to_line(const xy_t& x0, const xy_t& x1, const xy_t& x2)
{
    return std::abs((x2.x - x1.x) * (x1.y - x0.y) - (x1.x - x0.x) * (x2.y - x1.y))
           / x1.distance(x2);
}

static geo_t midpoint(const geo_t& p0, const geo_t& p1)
{
	/* Convert to Euclidean coordinates, compute the midpoint,
	 * convert back to spherical coordinates.
	 * This would work reasonably well if the points are not too
	 * far away from each other (i.e. leading to cancelling terms)
	 */
	double mx = 0.5 * (  std::cos(p0.lambda)*std::cos(p0.phi)
	                   + std::cos(p1.lambda)*std::cos(p1.phi));
	double my = 0.5 * (  std::sin(p0.lambda)*std::cos(p0.phi)
	                   + std::sin(p1.lambda)*std::cos(p1.phi));
	double mz = 0.5 * (std::sin(p0.phi) + std::sin(p1.phi));
	double nrm = std::sqrt(mx*mx + my*my + mz*mz);
	mx /= nrm;
	my /= nrm;
	mz /= nrm;
	mz = std::min(std::max(mz,-1.0), 1.0);
	return {std::atan2(my,mx), std::asin(mz)};
}


/*
 * This function refines, if required by the bisection criterion, a segment
 * ((x0,y0), (x1,y1)) until the minimum_node_distance boundary.
 * The bisection criterion ensures that round paths are approximated reasonably
 * well by the segments.
 *
 * The function furthermore tests whether any of the segments jump
 * a discontinuity. A discontinuity shows itself by nearly constant distance
 * in Euclidean space as the distance in geographic space gets smaller and
 * smaller.
 */


struct cut_t {
	xy_t first;
	xy_t second;

	cut_t(const xy_t& first, const xy_t& second) : first(first), second(second)
	{};
};

struct refined_segment_t {
	std::vector<path_xy_t> segments;
	std::vector<cut_t> discontinuities;
};

static xy_t spherical_stereographic(const geo_t& p, const geo_t& center,
                                    double a)
{
	/* Equations (21-2) to (21-4) from
	 * Snyder (1987) "Map Projections: A Working Manual".
	 * Assume k0 = 1.0 */
	const double sp1 = std::sin(center.phi);
	const double cp1 = std::cos(center.phi);
	const double sp = std::sin(p.phi);
	const double cp = std::cos(p.phi);
	const double cdl = std::cos(p.lambda - center.lambda);
	const double k = 2.0 / (1.0 + sp1 * sp + cp1 * cp * cdl);
	xy_t xy(a * k * cp * std::sin(p.lambda - center.lambda),
	        a * k * (cp1 * sp - sp1 * cp * cdl));
	return xy;
}


static refined_segment_t
refine_and_split_discontinuity(const std::vector<std::pair<geo_t,xy_t>>& path,
            const ProjWrapper& proj, double bisection_offset,
            double minimum_node_distance)
{
	/* The asymmetry constant.
	 * We compute (roughly) the geodesic center between two end points
	 * of a segment. Then we compute, in projected space, the Euclidean
	 * distances (d0,d1) between each of the two projected end points and
	 * the projected central points.
	 * If the smaller of the two, dmin, is smaller than the asymmetry constant
	 * times the larger of the two, CRITICAL_ASYMMETRY * dmax, we count this
	 * as a potential discontnuity. */
	constexpr double CRITICAL_ASYMMETRY = 0.1;

	constexpr double CA_SQ = CRITICAL_ASYMMETRY * CRITICAL_ASYMMETRY;

	auto potential_discontinuity = [=](double d2b0, double d2b1) -> bool {
		return std::min(d2b0,d2b1) < CA_SQ * std::max(d2b0,d2b1);
	};

	refined_segment_t refined;
	/* Stack of stacks.
	 * We start with one stack. Further stacks are added only if a
	 * discontinuity is detected and the path is split. */
	std::stack<std::stack<std::pair<geo_t,xy_t>>> remaining;
	remaining.emplace();
	for (auto it=path.rbegin(); it != path.rend(); ++it)
	{
		remaining.top().push(*it);
	}
	const double bisection_offset_2 = bisection_offset * bisection_offset;
	const double min_node_distance2 = minimum_node_distance
	                                  * minimum_node_distance;
	while (!remaining.empty()){
		std::pair<geo_t,xy_t> last(remaining.top().top());
		remaining.top().pop();
		refined.segments.emplace_back();
		refined.segments.back().push_back(last.second);
		while (!remaining.top().empty()){
			/* If the next node is closer than the minimum node distance, skip it
			 * (but only if its not the last node):
			 */
			const std::pair<geo_t,xy_t>& next = remaining.top().top();
			double distance2_non = last.second.distance_squared(next.second);
			if (distance2_non < min_node_distance2){
				if (remaining.size() == 1){
					/* Add the last node to the path: */
					refined.segments.back().push_back(next.second);
				}
				remaining.top().pop();
				continue;
			}

			/* Compute geographic center: */
			geo_t central(midpoint(last.first, next.first));
			xy_t xy_c(proj.project(central));

			/* We bisect only if the distance from the last to the central node
			 * is larger than the minimum distance (so as to prevent excess
			 * detail)
			 */
			double distance2_bisect = last.second.distance_squared(xy_c);

			/* Test whether we see asymmetry that might indicate a
			 * discontinuity: */
			double d2b_other = next.second.distance_squared(xy_c);
			if (potential_discontinuity(distance2_bisect,d2b_other)){
				/*
				 * Bisect the discontinuity.
				 */
				geo_t left, right;
				xy_t left_xy, right_xy;
				if (distance2_bisect > d2b_other){
					/* Discontinuity is between "last" and "central", that is,
					 * on the current path to central. */
					left = last.first;
					left_xy = last.second;
					right = central;
					right_xy = xy_c;
				} else {
					/* Discontinuity is between "central" and "next", that is,
					 * on the next path: */
					left = central;
					left_xy = xy_c;
					right = next.first;
					right_xy = next.second;
				}
				geo_t bisect_center(midpoint(left, right));
				xy_t stereo_l(spherical_stereographic(left, bisect_center,
				                                      proj.a()));
				xy_t stereo_r(spherical_stereographic(right, bisect_center,
				                                      proj.a()));
				while (stereo_l.distance_squared(stereo_r)
				       > 0.1 * min_node_distance2)
				{
					xy_t bice_xy(proj.project(bisect_center));
					double d2b_l = left_xy.distance_squared(bice_xy);
					double d2b_r = right_xy.distance_squared(bice_xy);
					if (d2b_l < CA_SQ * d2b_r){
						/* Need to bisect right interval. */
						left = bisect_center;
						left_xy = bice_xy;
					} else if (d2b_r < CA_SQ * d2b_l) {
						/* Need to bisect left interval. */
						right = bisect_center;
						right_xy = bice_xy;
					} else {
						/* No need to bisect! */
						break;
					}
					bisect_center = midpoint(left, right);
					stereo_l = spherical_stereographic(left, bisect_center,
					                                   proj.a());
					stereo_r = spherical_stereographic(right, bisect_center,
					                                   proj.a());
				}
				/* Now the line segment (left,right) is short enough for our
				 * minimum distance criterion and contains the discontinuity
				 * within.
				 * We split the remainder of the stack from "right" onwards
				 * to a new line segment / stack, and continue on the current
				 * stack (actually a new stack) until we reach "left". */
				remaining.top().push(std::make_pair(right, right_xy));

				/* Open the new stack: */
				remaining.emplace();
				remaining.top().push(std::make_pair(left, left_xy));

				/* Remember the discontinuity: */
				refined.discontinuities.emplace_back(left_xy, right_xy);
			} else {
				/*
				 * Test whether the bisection criterion is fulfilled.
				 * The bisection criterion is `bisection_offset`, the maximum
				 * offset that the central node can have from the line
				 * (last,next).
				 * This offset follows from trigonometry:
				 */
				if (   (distance2_non >= min_node_distance2)
					&& (distance2_bisect >= min_node_distance2)
					&& (distance_to_line(xy_c, next.second, last.second)
						>= bisection_offset_2))
				{
					/* Bisect, i.e. add the central point to the stack: */
					remaining.top().push(std::make_pair(central, xy_c));
				} else {
					/* Do not bisect, i.e. discard the central plot and add the
					 * stack top to the path: */
					refined.segments.back().push_back(next.second);
					last = next;
					remaining.top().pop();
				}
			}
		}
		remaining.pop();
	}

	return refined;
}



/*
 * This function cuts paths according to the following two criteria:
 *
 * (1) The path crosses the bounding box (xmin, xmax, ymin, ymax)
 *
 * (2) The path contains a discontinuity.
 *
 */
refined_path_t
flottekarte::crop_and_refine(const path_geo_t& geo_path,
    const ProjWrapper& proj, double xmin, double xmax, double ymin, double ymax,
    double bisection_offset, double minimum_node_distance)
{
	/* Helper methods: */
	auto xy_inside = [&](const xy_t& xy) -> bool {
		return (xy.x >= xmin) && (xy.x <= xmax) && (xy.y >= ymin)
		       && (xy.y <= ymax);
	};

	refined_path_t xy_paths;
	std::vector<std::pair<geo_t,xy_t>> path;
	size_t i=0;
	xy_t xy_i, xy_last;
	while (i<geo_path.size()){
		path.clear();

		/* First find the begin of a segment that is inside the map
		 * area: */
		while ((i < geo_path.size())
		       && (!xy_inside(xy_i = proj.project(geo_path[i]))))
		{
			++i;
			xy_last = xy_i;
		}
		if (i == geo_path.size())
			break;
		/* If we have found a point inside and i > 0, we have to find the
		 * point where the map boundary is intercepted. We do so later when
		 * we have refined the line according to the bisection
		 * criterion. For now, just additionally add the previous point.
		 */
		if (i > 0){
			path.emplace_back(geo_path[i-1].to_radians(), xy_last);
		}
		path.emplace_back(geo_path[i].to_radians(), xy_i);

		/* Find the end of the path: */
		while ((i < geo_path.size())
		       && (xy_inside(xy_i = proj.project(geo_path[i]))))
		{
			path.emplace_back(geo_path[i].to_radians(), xy_i);
			++i;
		}

		/* Add the last element (possibly), which is already outside the
		 * the map area. Again we will crop later after refinement: */
		if (i < geo_path.size()){
			path.emplace_back(geo_path[i].to_radians(), xy_i);
			xy_last = xy_i;
			++i;
		}

		/* Now we consider only the current entry in `paths`.
		 * First refine:
		 */
		refined_segment_t refined(
		    refine_and_split_discontinuity(path, proj, bisection_offset,
		                                   minimum_node_distance)
		);

		/* Add potentially existing cuts: */
		for (auto cut : refined.discontinuities){
			xy_paths.cuts.push_back(cut.first);
			xy_paths.cuts.push_back(cut.second);
		}

		/* Now consider all of the splits: */
		bool inside = false;
		std::optional<xy_t> sub_path_begin;
		std::vector<xy_t>::const_iterator it0, it1;
		for (const std::vector<xy_t>& path_refined : refined.segments){
			if (path_refined.size() == 0){
				std::cerr << "Encountered unexpectedly a zero-length path\n";
				throw std::runtime_error("Path problem.");
			}

			/* Now remove the possible outside nodes from the beginning and the
			 * end of the path. Also, the node bisection might have created
			 * nodes outside the map area. Hence, we traverse the `path_refined`
			 * and split it into paths [x(1), x(2), ..., x(n-1), x(n)], where
			 * x(1) and x(n) are outside (or the beginning/end of a
			 * discontinuity that is fully inside) and the remaining are inside.
			 */
			inside = xy_inside(path_refined[0]);
			it0 = path_refined.begin();
			for (auto it=path_refined.cbegin()+1; it != path_refined.cend();
			     ++it)
			{
				bool next_inside = xy_inside(*it);
				if (next_inside != inside){
					if (next_inside){
						/* Changing from outside to inside.
						 *
						 * Compute the intersection with the map boundary:
						 */
						sub_path_begin
						   = boundary_intersection(*(it-1), *it, xmin, xmax,
						                           ymin, ymax);

						/* Remember the first full inside node: */
						it0 = it;

						/* Remember the cut: */
						xy_paths.cuts.push_back(sub_path_begin.value());
					} else {
						/* Changing from inside to outside.
						 *
						 * Compute the intersection with the map boundary:
						 */
						xy_t sub_path_end(boundary_intersection(*it, *(it-1),
						                             xmin, xmax, ymin, ymax));

						/* Assemble the path: */
						xy_paths.segments.emplace_back();
						if (sub_path_begin)
							xy_paths.segments.back()
								    .push_back(sub_path_begin.value());
						xy_paths.segments.back().insert(
						    xy_paths.segments.back().cend(), it0, it
						);
						xy_paths.segments.back().push_back(sub_path_end);

						/* Remember the cut: */
						xy_paths.cuts.push_back(sub_path_end);

						/* Forget the sub path begin: */
						sub_path_begin.reset();
					}

					/* Set the new state. */
					inside = next_inside;
				}
			}

			/* Add the possible last path: */
			if (inside){
				/* Assemble the path: */
				xy_paths.segments.emplace_back();
				if (sub_path_begin)
					xy_paths.segments.back().push_back(sub_path_begin.value());
				xy_paths.segments.back().insert(xy_paths.segments.back().cend(),
				                                it0, path_refined.cend());
			}
		}
	}

	return xy_paths;
}
