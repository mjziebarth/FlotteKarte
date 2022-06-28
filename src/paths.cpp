/*
 * Path processing facilities.
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
using flottekarte::deg2rad;
using flottekarte::ProjWrapper;
using flottekarte::crop_and_refine;



static xy_t boundary_intersection(const xy_t& out, const xy_t& in, double xmin,
                                  double xmax, double ymin, double ymax)
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


static std::vector<xy_t>
refine_path(const std::vector<std::pair<geo_t,xy_t>>& path,
            const ProjWrapper& proj, double bisection_offset,
            double minimum_node_distance)
{
	std::vector<xy_t> path_refined;
	std::stack<std::pair<geo_t,xy_t>> remaining;
	for (auto it=path.rbegin(); it != path.rend(); ++it)
	{
		remaining.push(*it);
	}
	std::pair<geo_t,xy_t> last(remaining.top());
	remaining.pop();
	path_refined.push_back(last.second);
	const double bisection_offset_2 = bisection_offset * bisection_offset;
	const double min_node_distance2 = minimum_node_distance
	                                  * minimum_node_distance;
	while (!remaining.empty()){
		/* If the next node is closer than the minimum node distance, skip it
		 * (but only if its not the last node):
		 */
		const std::pair<geo_t,xy_t>& next = remaining.top();
		double distance2_non = last.second.distance_squared(next.second);
		if (distance2_non < min_node_distance2 && remaining.size() > 1){
			remaining.pop();
			continue;
		}

		/* Compute geographic center: */
		geo_t central(midpoint(last.first, next.first));
		xy_t xy_c(proj.project(central));

		/* We bisect only if the distance from the last to the central node
		 * is larger than the minimum distance (so as to prevent excess detail)
		 */
		double distance2_bisect = last.second.distance_squared(xy_c);

		/* Test whether the bisection criterion is fulfilled.
		 * The bisection criterion is `bisection_offset`, the maximum
		 * offset that the central node can have from the line (last,next).
		 * This offset follows from trigonometry:
		 */
		if (   (distance2_non >= min_node_distance2)
		    && (distance2_bisect >= min_node_distance2)
		    && (distance_to_line(xy_c, next.second, last.second)
		        >= bisection_offset_2))
		{
			/* Bisect, i.e. add the central point to the stack: */
			remaining.push(std::make_pair(central, xy_c));
		} else {
			/* Do not bisect, i.e. discard the central plot and add the
			 * stack top to the path: */
			path_refined.push_back(next.second);
			last = next;
			remaining.pop();
		}
	}

	return path_refined;
}


static double normed_vec_dot_product(const xy_t& ref, const xy_t& x0,
                                     const xy_t& x1)
{
	double vx0 = x0.x - ref.x;
	double vy0 = x0.y - ref.y;
	double vx1 = x1.x - ref.x;
	double vy1 = x1.y - ref.y;
	/* Product of the norms of the two vectors: */
	double nrm = std::sqrt((vx0*vx0 + vy0*vy0) * (vx1*vx1 + vy1*vy1));
	return (vx0 * vx1 + vy0*vy1) / nrm;
}


static void cut_at_angles(std::vector<path_xy_t>& paths,
                          const ProjWrapper& proj, double cut_at_angle_degrees,
                          bool cyclic)
{
	/* Testing criticality: */
	const double cos_cut_at_angle = std::cos(deg2rad(cut_at_angle_degrees));
	const bool try_cycle = cyclic && paths.size() == 1 && paths[0].size() > 2;

	auto is_critical = [&](const path_xy_t& path, size_t i) -> bool {
		size_t last, next;
		if (i == 0 || i == path.size()-1){
			if (try_cycle){
				if (i == 0){
					last = path.size()-1;
					next = 1;
				} else {
					last = path.size()-2;
					next = 0;
				}
			}
			return false;
		} else {
			last = i-1;
			next = i+1;
		}
		double dot_p = normed_vec_dot_product(path[i], path[last],
		                                      path[next]);
		return -dot_p < cos_cut_at_angle;
	};

	std::set<size_t> critical_nodes;
	std::vector<path_xy_t> append;
	for (path_xy_t& path : paths){
		const size_t M = path.size();
		if (M < 3){
			/* Cannot compute critiality for short paths: */
			continue;
		}
		critical_nodes.clear();

		/* Compute the mask that indicates whether the angle at each node
		 * is critial: */
		for (size_t i=1; i<M-1; ++i){
			if (is_critical(path, i))
				critical_nodes.insert(i);
		}

		/* Early continue if nothing is critical: */
		if (critical_nodes.size() == 0)
			continue;

		bool cycle = try_cycle;
		if (try_cycle){
			/* */
			bool end_critical = false;
			/* Now make sure that 0 != M-1 so as to be able to assess
			 * criticality: */
			if (path[0] == path.back()){
				path.pop_back();
			}
			if (is_critical(path, 0)){
				critical_nodes.insert(0);
				end_critical = true;
			}
			if (is_critical(path, path.size()-1)){
				critical_nodes.insert(path.size()-1);
				end_critical = true;
			}
			/* Cycle around only if the end is not critical: */
			cycle = !end_critical;

		} else {
			/* Handle the edge nodes. If possible, they take the blame if
			 * node 1 or M-1 are critical: */
			if (critical_nodes.contains(1) && !critical_nodes.contains(2)){
				/* Node 0 is the culprit: */
				critical_nodes.insert(0);
				critical_nodes.erase(1);
			}
			if (critical_nodes.contains(M-2) && !critical_nodes.contains(M-3)){
				/* Node M-1 is the culprit: */
				critical_nodes.insert(M-1);
				critical_nodes.erase(M-2);
			}
		}

		/* Check if there are enoug nodes remainig:: */
		const size_t remaining = path.size() - critical_nodes.size();
		if (remaining == 0){
			path.clear();
			continue;
		}

		/* Now create the reduced vector: */
		if (cycle){
			size_t begin = *critical_nodes.begin() + 1;
			for (auto it=++critical_nodes.cbegin(); it != critical_nodes.cend();
			     ++it)
			{
				append.emplace_back(path.cbegin()+begin, path.cbegin() + *it);
				begin = *it + 1;
			}
			if (begin < path.size() - 1)
				append.emplace_back(path.cbegin()+begin, path.cend());
			else
				append.emplace_back();

			/* Decide whether to insert the  */
			append.back().insert(append.back().cend(), path.cbegin(),
			                     path.cbegin() + *critical_nodes.begin());
			path.clear();
		} else {
			size_t begin = 0;
			for (size_t i : critical_nodes){
				append.emplace_back(path.cbegin()+begin, path.cbegin()+i);
				begin = i+1;
			}
			append.emplace_back(path.cbegin()+begin, path.cend());
			path.clear();
		}
	}

	/* Add the split paths: */
	paths.insert(paths.cend(), append.cbegin(), append.cend());

	/* Clear all empty paths: */
	auto to = paths.begin();
	for (auto from = paths.cbegin(); from != paths.cend(); ++from){
		if (from->size() > 1){
			*to = *from;
			++to;
		}
	}
	paths.resize(to - paths.begin());



	for (size_t j=0; j<paths.size(); ++j){
		const path_xy_t& p = paths[j];
		if (p.size() < 3)
			continue;
		for (size_t i=0; i<p.size(); ++i){
			if (is_critical(p, i)){
				std::cerr << "found a critical path after they should be "
				             "cleared!\n    i=" << i << " / "
				             << p.size() << "\n"
				             << "   try_cycle = " << try_cycle << "\n"
				             << "   j=" << j << " / " << paths.size() << "\n";
			}
		}
	}
}


/* This function aims to identify discontinuities.
 */
static void cut_discontinuities(std::vector<path_xy_t>& paths,
                                const ProjWrapper& proj)
{
	typedef std::vector<geo_t> geopath_t;
	std::stack<std::pair<path_xy_t,geopath_t>> split_off;

	/*  */
	auto find_discontinuity = [&](const path_xy_t& path,
	                              const geopath_t& gpath) -> size_t
	{
		for (size_t i=0; i<path.size()-1; ++i){
			/* Compute the center of the line segment (i, i+1) */
			geo_t center(midpoint(gpath[i],gpath[i+1]));
			xy_t c_xy(proj.project(center));

			/* Distances of both line segments: */
			double d0 = path[i].distance(c_xy);
			double d1 = path[i+1].distance(c_xy);

			/* Test for asymmetry: */
			if (d0 > 5 * d1 || d1 > 5 * d0){
				/* Split: */
				return i+1;
			}
		}
		return path.size();
	};

	for (path_xy_t& path : paths){
		if (path.size() < 2)
			continue;
		/* Obtain all the geographic coordinates: */
		std::vector<geo_t> geo_path(path.size());
		for (size_t i=0; i<path.size(); ++i){
			geo_path[i] = proj.inverse(path[i]);
		}

		/* Now iterate the lines, compute the centers in geographic coordinates,
		 * and see whether there is a crass asymmetry in the lengths of the
		 * two splits: */
		size_t disco = find_discontinuity(path, geo_path);
		if (disco != path.size()){
			/* Found a discontinuity. Remove it and add the second half to the
			 * stack: */
			std::cout << "found discontinuity at index " << disco << " / "
			          << path.size() << "\n" << std::flush;
			if (path.size() - disco > 1){
				path_xy_t tmp0(path.begin()+disco, path.end());
				geopath_t tmp1(geo_path.begin()+disco, geo_path.end());
				split_off.emplace(tmp0, tmp1);
			}
			path.resize(disco);
		}
	}

	/* Clear all empty paths: */
	auto to = paths.begin();
	for (auto from = paths.cbegin(); from != paths.cend(); ++from){
		if (from->size() > 1){
			*to = *from;
			++to;
		}
	}
	paths.resize(to - paths.begin());

	/* Now iterate all the split off discontinuities: */
	while (!split_off.empty()){
		const path_xy_t& path = split_off.top().first;
		const geopath_t& gpath = split_off.top().second;
		size_t disco = find_discontinuity(path, gpath);
		if (disco == path.size()){
			paths.push_back(path);
			split_off.pop();
		} else {
			/* Found a discontinuity. Remove it and add the second half to the
			 * stack: */
			if (disco > 1){
				paths.emplace_back(path.cbegin(), path.cbegin()+disco);
			}
			if (path.size() - disco > 1){
				path_xy_t tmp0(path.begin()+disco, path.end());
				geopath_t tmp1(gpath.begin()+disco, gpath.end());
				split_off.pop();
				split_off.emplace(tmp0, tmp1);
			} else {
				split_off.pop();
			}
		}
	}
}





std::vector<path_xy_t>
flottekarte::crop_and_refine(const path_geo_t& geo_path, const ProjWrapper& proj,
                          double xmin, double xmax, double ymin, double ymax,
                          double bisection_offset, double minimum_node_distance,
                          double cut_at_angle_degrees, bool cyclic)
{
	/* Helper methods: */
	auto xy_inside = [&](const xy_t& xy) -> bool {
		return (xy.x >= xmin) && (xy.x <= xmax) && (xy.y >= ymin)
		       && (xy.y <= ymax);
	};

	std::vector<path_xy_t> xy_paths;
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
		 * point where the map boundary is intercepted. We do so later,
		 * however, when we have refined the line according to the bisection
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
		std::vector<xy_t>
		path_refined(refine_path(path, proj, bisection_offset,
		                         minimum_node_distance));
		if (path_refined.size() == 0){
			std::cerr << "Encountered unexpectedly a zero-length path\n";
			throw std::runtime_error("Path problem.");
		}

		/* Now remove the possible outside nodes from the beginning and the
		 * end of the path. Also, the node bisection might have created
		 * nodes outside the map area. Hence, we traverse the `path_refined`
		 * and split it into paths [x(1), x(2), ..., x(n-1), x(n)], where
		 * x(1) and x(n) are outside and the remaining are inside.
		 */
		std::vector<xy_t>::const_iterator it0, it1;
		bool inside = xy_inside(path_refined[0]);
		it0 = path_refined.begin();
		std::optional<xy_t> sub_path_begin;
		for (auto it=path_refined.cbegin()+1; it != path_refined.cend(); ++it)
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
				} else {
					/* Changing from inside to outside.
					 *
					 * Compute the intersection with the map boundary:
					 */
					xy_t sub_path_end(boundary_intersection(*it, *(it-1),
					                             xmin, xmax, ymin, ymax));

					/* Assemble the path: */
					xy_paths.emplace_back();
					if (sub_path_begin)
						xy_paths.back().push_back(sub_path_begin.value());
					xy_paths.back().insert(xy_paths.back().cend(), it0, it);
					xy_paths.back().push_back(sub_path_end);

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
			xy_paths.emplace_back();
			if (sub_path_begin)
				xy_paths.back().push_back(sub_path_begin.value());
			xy_paths.back().insert(xy_paths.back().cend(), it0,
			                       path_refined.cend());
		}
	}

	/* Now handle the cutting criterion: */
	cut_at_angles(xy_paths, proj, cut_at_angle_degrees, cyclic);

	/* Finally, identify and cut discontinuities: */
	cut_discontinuities(xy_paths, proj);

	return xy_paths;
}
