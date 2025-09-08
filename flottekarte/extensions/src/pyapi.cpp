/*
 * API used for interfacing with Python.
 *
 * Authors: Malte J. Ziebarth (ziebarth@gfz-potsdam.de)
 *
 * Copyright (C) 2022 Deutsches GeoForschungsZentrum Potsdam,
 *                    Malte J. Ziebarth,
 *          2024-2025 Technische Universität München
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

#include <../include/pyapi.hpp>
#include <../include/invert.hpp>
#include <../include/projwrapper.hpp>
#include <../include/griddedinverter.hpp>
#include <../include/gradient.hpp>
#include <../include/tickfinder.hpp>
#include <../include/grid.hpp>
#include <../include/boundary.hpp>
#include <../include/streamlines.hpp>
#include <../include/azimuth.hpp>
#include <../include/objectbook.hpp>
#include <iostream>
#include <cmath>

using flottekarte::xy_t;
using flottekarte::geo_t;
using flottekarte::cut_t;
using flottekarte::axis_t;
using flottekarte::tick_t;
using flottekarte::modulo;
using flottekarte::rad2deg;
using flottekarte::deg2rad;
using flottekarte::path_xy_t;
using flottekarte::ProjError;
using flottekarte::geo_grid_t;
using flottekarte::ProjWrapper;
using flottekarte::geo_degrees_t;
using flottekarte::AugmentedProj;
using flottekarte::segment_tick_t;
using flottekarte::GriddedInverter;
using flottekarte::generate_grid_lines;
using flottekarte::gradient_descent_inverse_project;
using flottekarte::bounding_polygon;
using flottekarte::compute_ticks;
using flottekarte::Gradient;
using flottekarte::FORWARD_5POINT;
using flottekarte::ObjectBook;

/*void project_data(const char* proj_str, unsigned long Npoints,
                  const double* lon, const double* lat,
                  double* xy)
{

}*/

/*
 * Object types.
 */
struct grid_lines_t {
    std::vector<path_xy_t> paths;
    std::vector<cut_t> cuts;
};


/*
 * Caching the created objects:
 */
using ObjectBook_t = ObjectBook<
    size_t,
    AugmentedProj,
    grid_lines_t,
    path_xy_t,
    std::vector<path_xy_t>
>;

static ObjectBook_t object_book;



using index_t = uint64_t;





/*
 * Creating the ProjWrapper.
 */
int create_proj_wrapper_from_proj_str(
    const char* proj_str,
    index_t& oid
)
{
    try {
        /* Try to generate the PROJ wrapper: */
        oid = object_book.emplace(AugmentedProj(proj_str));
    } catch (const ProjError& err) {
        std::cout << "error occurred.\n what: " << err.what() << "\n"
                  << std::flush;
        return 1;
    } catch (...) {
        #ifdef DEBUG
        std::cout << "error occurred.\n" << std::flush;
        #endif
        return 2;
    }
    return 0;
}




int inverse_project_data_optimize(
    index_t oid,
    unsigned long Npoints,
    const double* x,
    const double* y,
    double* lon_lat
)
{
    /* Initialize the projection: */
    #ifdef DEBUG
    std::cout << "initializing proj wrapper\n" << std::flush;
    #endif
    try {
        std::shared_ptr<AugmentedProj> pptr(object_book.get<AugmentedProj>(oid));
        ProjWrapper& proj(*pptr);
        GriddedInverter ginv(proj, 100, 50);

        /* Parallel projection: */
        #ifdef DEBUG
        std::cout << "start projection loop...\n" << std::flush;
        #endif

        #pragma omp parallel for
        for (size_t i=0; i<Npoints; ++i){
            xy_t xy(x[i], y[i]);
            try {
                /* First try to use PROJ inverse of projection: */
                geo_t lola(proj.inverse(xy));
                lon_lat[2*i] = rad2deg(lola.lambda);
                lon_lat[2*i+1] = rad2deg(lola.phi);
            } catch (const ProjError& e) {
                /* If this fails, invert with gradient descent: */
                geo_t start(ginv(xy));
                geo_t lola(gradient_descent_inverse_project(proj, xy,
                                                            start.lambda,
                                                            start.phi));
                lon_lat[2*i] = rad2deg(lola.lambda);
                lon_lat[2*i+1] = rad2deg(lola.phi);
            }
        }

        #ifdef DEBUG
        std::cout << "cleaning up!\n" << std::flush;
        #endif

        /* Success: */
        return 0;
    } catch (const ProjError& err) {
        std::cout << "error occurred.\n what: " << err.what() << "\n"
                  << std::flush;
        return 1;
    } catch (...) {
        #ifdef DEBUG
        std::cout << "error occurred.\n" << std::flush;
        #endif
        return 1;
    }
}


int gradients_east_north(
    index_t oid,
    unsigned long Npoints,
    const double* lon,
    const double* lat,
    double* gradient_east,
    double* gradient_north,
    double stencil_delta
)
{
    try {
        /* Initialize the projection: */
        std::shared_ptr<AugmentedProj> pptr(object_book.get<AugmentedProj>(oid));
        ProjWrapper& proj(*pptr);

        /* Parallel evaluation of the gradients: */
        #pragma omp parallel for
        for (size_t i=0; i<Npoints; ++i){
            /* Compute coordinate gradients in east and north: */
            geo_t lola({deg2rad(lon[i]), deg2rad(lat[i])});
            Gradient<FORWARD_5POINT> gradient(proj, lola, stencil_delta);

            /* Save the values: */
            gradient_east[2*i]    = gradient.gx_east();
            gradient_east[2*i+1]  = gradient.gy_east();
            gradient_north[2*i]   = gradient.gx_north();
            gradient_north[2*i+1] = gradient.gy_north();
        }

        /* Success: */
        return 0;
    } catch (...) {
        return 1;
    }
}


int scale_factors(
    index_t oid,
    unsigned long Npoints,
    const double* lon,
    const double* lat,
    double* kh,
    double stencil_delta
)
{
    try {
        /* Initialize the projection: */
        std::shared_ptr<AugmentedProj> pptr(object_book.get<AugmentedProj>(oid));
        ProjWrapper& proj(*pptr);

        /* Parallel evaluation of the gradients: */
        #pragma omp parallel for
        for (size_t i=0; i<Npoints; ++i){
            /* Compute coordinate gradients in east and north: */
            geo_t lola({deg2rad(lon[i]), deg2rad(lat[i])});
            Gradient<FORWARD_5POINT> g(proj, lola, stencil_delta);

            /* Save the values: */
            const double k = std::sqrt(  g.gx_east() * g.gx_east()
                                       + g.gy_east() * g.gy_east());
            kh[2*i]   = k;
            const double h = std::sqrt(  g.gx_north() * g.gx_north()
                                       + g.gy_north() * g.gy_north());
            kh[2*i+1] = h;
        }

        /* Success: */
        return 0;
    } catch (...) {
        return 1;
    }
}


int compute_axes_ticks(
    index_t oid,
    size_t Nseg,
    const double* vertices,
    double tick_spacing_degree,
    unsigned int max_ticks,
    unsigned int* segments,
    double* tick_vertices,
    unsigned char* which_ticks,
    unsigned int* Nticks
)
{
    /* Initialize the projection: */
    #ifdef DEBUG
    std::cout << "initializing proj wrapper\n" << std::flush;
    #endif
    try {
        std::shared_ptr<AugmentedProj> pptr(object_book.get<AugmentedProj>(oid));
        AugmentedProj& proj(*pptr);

        /* Parallel projection: */
        #ifdef DEBUG
        std::cout << "tick creation loop...\n" << std::flush;
        #endif

        /* Prepare everything for loop execution: */
        constexpr tick_t tick[2] = {flottekarte::TICK_LON,
                                    flottekarte::TICK_LAT};
        path_xy_t path(Nseg);
        for (xy_t& xy : path){
            xy.x = *vertices;
            ++vertices;
            xy.y = *vertices;
            ++vertices;
        }

        const size_t Nmax = static_cast<size_t>(max_ticks);
        std::array<std::vector<segment_tick_t>,2> ticks;
        for (int i=0; i<2; ++i){
            ticks[i]
               = compute_ticks(proj, tick[i], path, tick_spacing_degree);
        }
        size_t n=0;
        for (int i=0; i<2; ++i){
            const size_t Nremain = Nmax - n;
            const size_t Ni = std::min<size_t>(ticks[i].size(), Nremain);
            for (size_t j=0; j<Ni; ++j){
                *tick_vertices = ticks[i][j].tick.lon;
                ++tick_vertices;
                *tick_vertices = ticks[i][j].tick.lat;
                ++tick_vertices;
                *segments = ticks[i][j].segment;
                ++segments;
                *which_ticks = static_cast<unsigned char>(tick[i]);
                ++which_ticks;
            }
            n += Ni;
            if (n == Nmax){
                std::cout << "Warning: Reached maximum number of ticks.\n";
                break;
            }
        }
        *Nticks = static_cast<unsigned int>(n);

        #ifdef DEBUG
        std::cout << "cleaning up!\n" << std::flush;
        #endif

        /* Success: */
        return 0;

    } catch (const ProjError& err) {
        std::cout << "error occurred.\n what: " << err.what() << "\n"
                  << std::flush;
        return 1;
    } catch (...) {
        #ifdef DEBUG
        std::cout << "error occurred.\n" << std::flush;
        #endif
        return 1;
    }
}


/******************************************************************************
 *                             Azimuth conversion                             *
 ******************************************************************************/

int azimuth_geographic_to_local_on_grid_inplace(
    index_t oid,
    double xmin,
    double xmax,
    size_t nx,
    double ymin,
    double ymax,
    size_t ny,
    double* azimuth_rad,
    size_t Nazi,
    double stencil_delta
)
{
    try {
        /* Initialize the projection: */
        std::shared_ptr<AugmentedProj> pptr(object_book.get<AugmentedProj>(oid));
        ProjWrapper& proj(*pptr);

        /* Grid description: */
        const long double dx = (xmax - (long double)xmin) / nx;
        const long double dy = (ymax - (long double)ymin) / ny;

        /* Parallel evaluation of the gradients: */
        #pragma omp parallel for
        for (size_t k=0; k<Nazi; ++k){
            /* Compute the geographic coordinates of the grid point: */
            size_t i = k / ny;
            size_t j = k % ny;
            xy_t xy(xmin + i*dx, ymin + j * dy);

            /* Compute coordinate gradients in east and north.
             * Note that 'Gradient' computes the derivative of the projected
             * coordinates by true physical space, i.e. along two equal
             * infinitesimal vectors pointing in the north (lat) and east (lon)
             * direction.
             * Hence, the gradient quantifies the local affine transform of
             * true rectangular coordinates. We can hence add the gradients
             * scaled by the angular component to retrieve an unnormed
             * directional vector in projected space.
             */
            geo_t lola(proj.inverse(xy));
            Gradient<FORWARD_5POINT> g(proj, lola, stencil_delta);

            /* Directional vectors: */
            xy_t east(g.gx_east(), g.gy_east());
            xy_t north(g.gx_north(), g.gy_north());

            /* Now find the azimuth direction: */
            xy_t azi_dir(
                std::sin(azimuth_rad[k]) * east
                + std::cos(azimuth_rad[k]) * north
            );

            /* Find the local azimuth of that vector: */
            azimuth_rad[k] = std::atan2(azi_dir.x, azi_dir.y);
        }

        /* If all Nazi are converted, we're done! */
        return 0;
    } catch (...) {
        return 1;
    }
}


int unwrap_azimuth_field(
    double* angle,
    uint32_t nx,
    uint32_t ny,
    size_t Nmax,
    double cost_beta
)
{
    try {
        flottekarte::unwrap_azimuth_field(
            angle, nx, ny, Nmax, cost_beta
        );
    } catch (const std::exception& e){
        std::cerr << "Error in unwrap_and_smooth_azimuth_field: "
                  << e.what() << "\n";
        return 1;
    } catch (...) {
        return 2;
    }

    return 0;
}


/******************************************************************************
 *                                 Grid lines                                 *
 ******************************************************************************/


/*
 * This is the first of a two-part function. Computes grid lines according to
 * settings (projection, x- & ylim, tick spacing, bisection offset,
 * minimum distance between path-adjacent nodes) and returns 1) the length
 * of the resulting path (Npath) and 2) a pointer to the structure holding
 * all the required information (struct_ptr).
 * It has to be followed, after allocating two suitable numpy arrays,
 * by a call to save_grid_lines, which transfers the path contained in
 * struct_ptr and frees the allocated space. If this call is omitted, a memory
 * leak will occur. If this call a second time, unallocated memory will be
 * accessed.
 */
int compute_grid_lines(
    index_t proj_oid,
    double xmin,
    double xmax,
    double ymin,
    double ymax,
    int tick_spacing_degree,
    double bisection_offset,
    double minimum_node_distance,
    double max_lat,
    index_t* gridlines_oid,
    size_t* Npath,
    size_t* Ncut
)
{
    /* Empty initialization: */
    if (!gridlines_oid){
        std::cerr << "No pointer to write gridlines struct object ID given.\n";
        return 3;
    }
    if (!Npath){
        std::cerr << "No pointer to write Npath given.\n";
        return 4;
    }
    if (!Ncut){
        std::cerr << "No pointer to write Ncut given.\n";
        return 5;
    }
    *gridlines_oid = 0;
    *Npath = 0;
    *Ncut = 0;

    try {
        /* Create the projection wrapper: */
        std::shared_ptr<AugmentedProj> pptr(object_book.get<AugmentedProj>(proj_oid));
        ProjWrapper& proj(*pptr);

        /* Generate the grid lines structure: */
        grid_lines_t glines;

        /* Compute the grid lines: */
        geo_grid_t grid(
            generate_grid_lines(
                proj,
                xmin,
                xmax,
                ymin,
                ymax,
                tick_spacing_degree,
                bisection_offset,
                minimum_node_distance,
                max_lat
            )
        );
        glines.paths.swap(grid.paths);
        glines.cuts.swap(grid.cuts);

        /* Compute the number of grid lines: */
        size_t npath = 0;
        for (const path_xy_t& path : glines.paths){
            npath += path.size();
        }
        size_t ncuts = glines.cuts.size();

        /* Move the grid lines struct: */
        *gridlines_oid = object_book.emplace<grid_lines_t>(std::move(glines));

        /* Save the sizes:: */
        *Npath = npath;
        *Ncut = ncuts;

        /* Success. */
        return 0;

    } catch (const ProjError& err){
        /* Could not create the projection object: */
        std::cerr << "ProjError: " << err.what() << "\n";
        return 2;
    } catch (const std::bad_alloc& err){
        /* Possibly could not allocate the grid lines object. */
        std::cerr << "Memory allocation error.\n";
        return 6;
    }
}

constexpr uint8_t MPL_MOVETO = 1;
constexpr uint8_t MPL_LINETO = 2;

/*
 * Second part of a two-part function. Call exactly one time after
 * compute_grid_lines has successfully completed.
 */
int save_grid_lines(
    index_t oid,
    double* vertices,
    uint8_t* codes,
    double* cut_vertices,
    uint8_t* cut_axes,
    double* cut_coords
)
{
    try {
        /* Get the grid lines object: */
        std::shared_ptr<grid_lines_t> pptr(object_book.get<grid_lines_t>(oid));
        grid_lines_t& glines(*pptr);

        /* Fill the vertices and codes arrays: */
        for (const path_xy_t& path : glines.paths){
            for (size_t i=0; i<path.size(); ++i){
                /* Set the code: */
                if (i == 0){
                    *codes = MPL_MOVETO;
                } else {
                    *codes = MPL_LINETO;
                }

                /* Set the vertices: */
                *vertices = path[i].x;
                ++vertices;
                *vertices = path[i].y;

                /* Advance: */
                ++codes;
                ++vertices;
            }
        }

        /* Fill the cut points: */
        for (const cut_t& cut : glines.cuts){
            *cut_vertices = cut.point.x;
            ++cut_vertices;
            *cut_vertices = cut.point.y;
            ++cut_vertices;
        }
        for (const cut_t& cut : glines.cuts){
            *cut_axes = static_cast<uint8_t>(cut.tick_type);
            ++cut_axes;
        }
        for (const cut_t& cut : glines.cuts){
            *cut_coords = cut.coordinate;
            ++cut_coords;
        }

        /* Success. */
        return 0;
    } catch (const std::exception& err){
        std::cerr << "Error in save_gridlines: '"
            << err.what() << "'.\n" << std::flush;
        return 1;
    } catch (...){
        std::cerr << "Unknown error in save_gridlines.\n" << std::flush;
        return 2;
    }
}


int clean_grid_lines_struct(
    index_t oid
)
{
    try {
        object_book.erase(oid);

        return 0;
    } catch (const std::exception& err){
        std::cerr << "Error in clean_grid_lines_struct: '"
            << err.what() << "'.\n" << std::flush;
        return 1;
    } catch (...){
        std::cerr << "Unknown error in clean_grid_lines_struct.\n" << std::flush;
        return 2;
    }
}


/******************************************************************************
 *                   Computing the bounding polygon.                          *
 ******************************************************************************/

int compute_bounding_polygon(
    index_t oid,
    double xmin,
    double xmax,
    double ymin,
    double ymax,
    double atol,
    double bisection_offset,
    double minimum_node_distance,
    index_t* polygon_oid_ptr,
    size_t* Nvert
)
{
    /* Empty initialization: */
    if (!polygon_oid_ptr){
        std::cerr << "No pointer to write struct given.\n";
        return 3;
    }
    if (!Nvert){
        std::cerr << "No pointer to write Nver given.\n";
        return 4;
    }
    *Nvert = 0;

    try {
        /* Create the projection wrapper: */
        std::shared_ptr<AugmentedProj> pptr(object_book.get<AugmentedProj>(oid));
        AugmentedProj& proj(*pptr);

        /* Generate the bounding polygon structure: */
        path_xy_t poly;

        /* Compute the polygon: */
        bounding_polygon(
            proj,
            xmin,
            xmax,
            ymin,
            ymax,
            atol,
            bisection_offset,
            minimum_node_distance,
            poly
        );

        /* Number of vertices: */
        *Nvert = poly.size();

        *polygon_oid_ptr = object_book.emplace(std::move(poly));


        /* Success. */
        return 0;

    } catch (const ProjError& err){
        /* Could not create the projection object: */
        std::cerr << "ProjError: " << err.what() << "\n";
        return 2;
    } catch (const std::bad_alloc& err){
        /* Could not allocate memory somewhere. */
        std::cerr << "Memory allocation error in compute_bounding_polygon.\n";
        return 5;
    } catch (const std::exception& err){
        std::cerr << "Exception in compute_bounding_polygon: '"
            << err.what() <<  "'\n";
        return 6;
    } catch (...){
        std::cerr << "Unknown error in compute_bounding_polygon.\n" << std::flush;
        return 7;
    }
}


int save_bounding_polygon(
    index_t oid,
    double* vertices,
    double* angles
)
{
    try {
        std::shared_ptr<path_xy_t> pptr(object_book.get<path_xy_t>(oid));
        const path_xy_t& poly = *pptr;

        /* Check for empty polygon: */
        if (poly.empty())
            return 0;

        /* Sanity check: */
        if (!vertices){
            std::cerr << "No pointer to vertices given.\n";
            return 2;
        }
        if (!angles){
            std::cerr << "No pointer to angles given.\n";
            return 2;
        }

        /* Computing a segment angle: */
        auto segment_angle = [](const xy_t& last, const xy_t& next) -> double {
            /* First the forward direction this segment */
            const double dx = next.x - last.x;
            const double dy = next.y - last.y;
            /* The orthogonal outwards-pointing direction is
            * dx' = dy
            * dy' = -dx
            * This makes use of the fact that the bounding polygon winds
            * counterclockwise. Compute the angle of this direction in the
            * coordinate system used in the Python code: clockwise from the
            * positive y-axis.
            */
            const double angle = -modulo(rad2deg(std::atan2(-dx, dy)) - 90.0,
                                        360.0);
            if (angle < -180.0)
                return angle + 360.0;
            return angle;
        };

        /* Fill the arrays: */
        auto it = poly.cbegin();
        *vertices = it->x;
        ++vertices;
        *vertices = it->y;
        ++vertices;
        auto last = it;
        for (++it; it != poly.cend(); ++it){
            /* Fill the vertex: */
            *vertices = it->x;
            ++vertices;
            *vertices = it->y;
            ++vertices;
            /* Compute and fill the angle. */
            *angles = segment_angle(*last, *it);
            ++angles;
            last = it;
        }
        /* The last segment back to the start: */
        *angles = segment_angle(*last, poly[0]);

        return 0;
    } catch (const std::exception& err){
        std::cerr << "Error in save_bounding_polygon: '"
            << err.what() << "'.\n" << std::flush;
        return 1;
    } catch (...){
        std::cerr << "Unknown error in save_bounding_polygon.\n"
            << std::flush;
        return 2;
    }
}


int clean_bounding_polygon_struct(
    void* struct_ptr
)
{
    if (struct_ptr == nullptr)
        return 1;

    /* Delete the grid_lines_t struct, calling all necessary destructors: */
    path_xy_t* poly = static_cast<path_xy_t*>(struct_ptr);
    delete poly;

    /* Success. */
    return 0;
}

/******************************************************************************
 *                   Computing the streamline polygons.                       *
 ******************************************************************************/

/* Save all streamline polygons in a map. */


int compute_streamlines(
    double xmin,
    double xmax,
    size_t nx,
    double ymin,
    double ymax,
    size_t ny,
    const double* z,
    size_t Nz,
    double r,
    double ds_min,
    double width_scale,
    double epsilon,
    uint8_t width_mode,
    size_t* struct_id
)
{
    auto nan_or_inf = [](double x) -> bool {
        return std::isnan(x) || std::isinf(x);
    };

    /* Input sanity: */
    if (Nz != 3 * nx * ny)
        return 1;
    if (!z)
        return 2;
    if (!struct_id)
        return 3;
    if (nan_or_inf(xmin) || nan_or_inf(xmax) || nan_or_inf(ymin)
        || nan_or_inf(ymax))
        return 4;
    if (xmin >= xmax || ymin >= ymax)
        return 5;
    if (r >= 0.25 * std::min(xmax-xmin, ymax - ymin))
        // We will not actually have a visually pleasing amount of streamlines.
        return 6;
    if (nx < 2 || ny < 2)
        return 7;

    try {
        size_t oid = object_book.emplace(
            flottekarte::streamlines(
                xmin, xmax, nx, ymin, ymax, ny, z, r, ds_min, width_scale,
                epsilon, static_cast<flottekarte::width_mode_t>(width_mode)
            )
        );

        /* Success, save the object id: */
        *struct_id = oid;
    } catch (const std::exception& e){
        std::cerr << "exception: '" << e.what() << "'.\n" << std::flush;
        return 8;
    } catch (...) {
        return 8;
    }

    return 0;
}


size_t get_streamline_polygon_count(
    size_t oid
)
{
    try {
        return object_book.get<std::vector<path_xy_t>>(oid)->size();
    } catch (const std::exception& e){
        std::cerr << "exception: '" << e.what() << "'\n" << std::flush;
    } catch (...){
        std::cerr << "Unknown exception.\n" << std::flush;
    }
    return 0;
}


size_t get_streamline_polygon_size(
    size_t oid,
    size_t poly_id
)
{
    try {
        return object_book.get<std::vector<path_xy_t>>(oid)->at(poly_id).size();
    } catch (const std::exception& e){
        std::cerr << "exception: '" << e.what() << "'\n" << std::flush;
    } catch (...){
        std::cerr << "Unknown exception.\n" << std::flush;
    }
    return 0;
}


int save_streamline_polygon(
    size_t oid,
    size_t poly_id,
    double* out_xy,
    size_t Nout
)
{
    if (!out_xy)
        return 1;

    try {
        /* Get the correct path: */
        std::shared_ptr<std::vector<path_xy_t>> pptr
            = object_book.get<std::vector<path_xy_t>>(oid);
        const path_xy_t& path(pptr->at(poly_id));

        /* Now we found the correct polygon. Ensure that the buffer suffices: */
        if (2*path.size() != Nout)
            return 4;

        /* Now output the coordinates: */
        for (const xy_t& xy : path){
            *out_xy = xy.x;
            ++out_xy;
            *out_xy = xy.y;
            ++out_xy;
        }

        return 0;
    } catch (const std::exception& e){
        std::cerr << "exception: '" << e.what() << "'\n" << std::flush;
    } catch (...){
        std::cerr << "Unknown exception.\n" << std::flush;
    }

    return 1;
}


int delete_streamline_struct(
    size_t oid
)
{
    try {
        object_book.erase(oid);

        return 0;
    } catch (const std::exception& e){
        std::cerr << "exception: '" << e.what() << "'\n" << std::flush;
    } catch (...){
        std::cerr << "Unknown exception.\n" << std::flush;
    }

    return 1;
}