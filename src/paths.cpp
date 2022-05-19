


#include <../include/paths.hpp>



#include <optional>
#include <stack>
#include <iostream>
#include <cmath>
#include <set>


using projplot::xy_t;
using projplot::geo_t;
using projplot::geo_degrees_t;
using projplot::path_xy_t;
using projplot::deg2rad;
using projplot::ProjWrapper;
using projplot::crop_and_refine;



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
		/* If the next node is closer than the minimum node distance, skip it:
		 */
		const std::pair<geo_t,xy_t>& next = remaining.top();
		double distance2_non = last.second.distance_squared(next.second);
		if (distance2_non < min_node_distance2){
			remaining.pop();
			continue;
		}

		/* Compute geographic center: */
		double lam_c = 0.5 * (last.first.lambda + next.first.lambda);
		double phi_c = 0.5 * (last.first.phi + next.first.phi);
		geo_t central({lam_c, phi_c});
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


std::vector<path_xy_t>
projplot::crop_and_refine(const path_geo_t& geo_path, const ProjWrapper& proj,
                          double xmin, double xmax, double ymin, double ymax,
                          double bisection_offset, double minimum_node_distance)
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

	return xy_paths;
}
