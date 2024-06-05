/*
 * Grid index.
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

#ifndef FLOTTEKARTE_INDEX_HPP
#define FLOTTEKARTE_INDEX_HPP

#include <stdexcept>

namespace flottekarte {

template<typename Integer>
class GridIndex
{
public:
    Integer i;
    Integer j;
    Integer ny;

    GridIndex(Integer i, Integer j, Integer nx, Integer ny)
       : i(i), j(j), ny(ny)
    {
        if (i+1 > nx || j+1 > ny)
            throw std::runtime_error("Grid index out of bounds.");
    }

    constexpr GridIndex() noexcept : i(0), j(0), ny(0)
    {}

    Integer flat() const
    {
        return ny * i + j;
    }

    GridIndex<Integer> next_y() const
    {
        return GridIndex<Integer>(i, j+1, ny);
    }

    GridIndex<Integer> previous_y() const
    {
        return GridIndex<Integer>(i, j-1, ny);
    }

    GridIndex<Integer> next_x() const
    {
        return GridIndex<Integer>(i+1, j, ny);
    }

    GridIndex<Integer> previous_x() const
    {
        return GridIndex<Integer>(i-1, j, ny);
    }

    GridIndex<Integer> next_xy() const
    {
        return GridIndex<Integer>(i+1, j+1, ny);
    }

    static GridIndex<Integer> from_linear(size_t k, Integer nx, Integer ny)
    {
        return GridIndex<Integer>(k / ny, k % ny, ny);
    }

protected:
    GridIndex(Integer i, Integer j, Integer ny) : i(i), j(j), ny(ny)
    {}
};



} // end namespace


#endif