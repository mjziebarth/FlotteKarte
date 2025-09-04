/*
 * Two-dimensional interpolation.
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

#include <array>
#include <limits>
#include <cmath>
#include <../include/types.hpp>
#include <../include/index.hpp>

#ifndef FLOTTEKARTE_INTERPOLATE_HPP
#define FLOTTEKARTE_INTERPOLATE_HPP

namespace flottekarte {

/*
 * Structure to represent a data grid
 * ----------------------------------
 */
template<bool copy_z_buffer>
struct grid2d_t;


template<>
struct grid2d_t<true>
{
public:
    double xmin;
    double xmax;
    double ymin;
    double ymax;
    size_t nx;
    size_t ny;

    grid2d_t(
        double xmin,
        double xmax,
        size_t nx,
        double ymin,
        double ymax,
        size_t ny,
        const double* z
    );

private:
    std::vector<double> data;

public:
    const double* z;
};


template<>
struct grid2d_t<false>
{
public:
    double xmin;
    double xmax;
    double ymin;
    double ymax;
    size_t nx;
    size_t ny;

    grid2d_t(
        double xmin,
        double xmax,
        size_t nx,
        double ymin,
        double ymax,
        size_t ny,
        const double* z
    );

public:
    const double* z;
};


/*
 * A tiling of the two-dimensional space.
 * --------------------------------------
 */
template<typename Integer>
class TileIndex : public GridIndex<Integer>
{
public:
    TileIndex(
        const xy_t& xy,
        double xmin,
        double xmax,
        Integer nx,
        double ymin,
        double ymax,
        Integer ny
    )
      : GridIndex<Integer>(
            std::min<Integer>(std::max<Integer>(
                    std::floor((nx-1) * (xy.x - xmin) / (xmax - xmin)),
                    0
                ),
                nx-1
            ),
            std::min<Integer>(std::max<Integer>(
                    std::floor((ny-1) * (xy.y - ymin) / (ymax - ymin)),
                    0
                ),
                ny-1
            ),
            ny
        )
    {
        if (xy.x > xmax || xy.y > ymax || xy.x < xmin || xy.y < ymin)
            throw std::runtime_error("Coordinate out of bounds in TileIndex.");
        if (GridIndex<Integer>::i + 1 > nx || GridIndex<Integer>::j+1 > ny)
            throw std::runtime_error("TileIndex corrupted.");
    }

    TileIndex(GridIndex<Integer>&& gid) : GridIndex<Integer>(gid)
    {}

    TileIndex(const GridIndex<Integer>& gid) : GridIndex<Integer>(gid)
    {}
};

/*
 * A two-dimensional linear interpolator of gridded data.
 * ------------------------------------------------------
 */
template<bool copy_z_buffer, uint_fast8_t D=1>
class LinearGrid2DInterpolator
{
public:
    LinearGrid2DInterpolator(
        double xmin,
        double xmax,
        size_t nx,
        double ymin,
        double ymax,
        size_t ny,
        const double* z
    ) : grid(xmin, xmax, nx, ymin, ymax, ny, z)
    {
        if (nx == 0 || ny == 0)
            throw std::runtime_error("'nx' and 'ny' need to be positive.");
    }


    std::array<double,D> operator()(const xy_t& xy) const
    {
        /* Sanity: */
        if (xy.x < grid.xmin || xy.x > grid.xmax || xy.y < grid.ymin
            || xy.y > grid.ymax)
        {
            std::array<double,D> res;
            for (uint_fast8_t i=0; i<D; ++i)
                res[i] = std::numeric_limits<double>::quiet_NaN();
            return res;
        }

        /* Get the index of the tile we're in: */
        TileIndex id(xy, grid.xmin, grid.xmax, grid.nx, grid.ymin,
                     grid.ymax, grid.ny);

        /* Get the two to four data: */
        const long double dx = (grid.xmax - static_cast<long double>(grid.xmin))
                                / (grid.nx-1);
        const long double dy = (grid.ymax - static_cast<long double>(grid.ymin))
                                / (grid.ny-1);
        std::array<double,D> result;
        const size_t i0 = D*id.flat();
        if (id.i == grid.nx - 1){
            if (id.j == grid.ny - 1){
                /* xy is exactly the upper right boundary: */
                std::copy(grid.z + i0, grid.z + i0 + D, result.begin());
            } else {
                /* xy is on the right boundary. */
                long double wy = (xy.y - (grid.ymin + id.j * dy)) / dy;

                const size_t iny = D*id.next_y().flat();
                for (uint_fast8_t i=0; i<D; ++i)
                    result[i] =  (1.0 - wy) * grid.z[i0+i]
                                + wy * grid.z[iny+i];
            }
        }
        else if (id.j == grid.ny - 1)
        {
            /* xy is on the top boundary. */
            long double wx = (xy.x - (grid.xmin + id.i * dx)) / dx;

            const size_t inx = D*id.next_x().flat();
            for (uint_fast8_t i=0; i<D; ++i)
                result[i] =  (1.0 - wx) * grid.z[i0+i]
                            + wx * grid.z[inx+i];
        }
        else
        {
            /* The general case: xy is not on the boundary. */
            long double wx = (xy.x - (grid.xmin + id.i * dx)) / dx;
            long double wy = (xy.y - (grid.ymin + id.j * dy)) / dy;

            const size_t iny = D*id.next_y().flat();
            const size_t inx = D*id.next_x().flat();
            const size_t inxy = D*id.next_xy().flat();

            for (uint_fast8_t i=0; i<D; ++i)
                result[i] = (1.0 - wx) * (1.0 - wy) * grid.z[i0+i]
                    + (1.0 - wx) * wy * grid.z[iny+i]
                    + wx * (1.0 - wy) * grid.z[inx+i]
                    + wx * wy * grid.z[inxy+i];
        }

        return result;
    }

private:
    grid2d_t<copy_z_buffer> grid;
};


}

#endif