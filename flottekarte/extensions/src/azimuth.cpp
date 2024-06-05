/*
 * A method to disentangle an azimuth field, hopefully yielding a pleasantly
 * smoothed field.
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

#include <../include/azimuth.hpp>
#include <../include/index.hpp>
#include <../include/types.hpp>


#include <cmath>
#include <array>
#include <vector>
#include <queue>
#include <limits>
#include <numbers>
#include <stddef.h>

namespace flottekarte {

typedef GridIndex<uint32_t> index_t;


static xy_t global_mean_velocity(
    const double* angle,
    const std::vector<int8_t>& sign,
    const size_t nx, const size_t ny)
{
    /*
     * Computes the global mean velocity.
     */
    long double vgx = 0.0;
    long double vgy = 0.0;
    for (size_t k=0; k<nx*ny; ++k){
        long double angle_ij = angle[k];
        int8_t s_ij = sign[k];
        vgx += s_ij * std::sin(angle_ij);
        vgy += s_ij * std::cos(angle_ij);
    }

    /* Normalize: */
    long double vgn = std::sqrt(vgx*vgx + vgy*vgy);

    return xy_t(vgx / vgn, vgy / vgn);
}


static double cost_contribution(
    GridIndex<uint32_t> ij,
    const double* angle,
    const std::vector<int8_t>& sign,
    const xy_t& v_global,
    uint32_t nx,
    double cost_beta
)
{
    /*
     * Computes the contribution of a single node to the global cost.
     */
    const uint32_t i = ij.i;
    const uint32_t j = ij.j;

    /* Check which neighbors we have: */
    std::array<GridIndex<uint32_t>,4> neighbors;
    uint_fast8_t Nnb = 0;
    if (ij.i > 0){
        neighbors[0] = ij.previous_x();
        ++Nnb;
    }
    if (i < nx-1){
        neighbors[Nnb] = ij.next_x();
        ++Nnb;
    }
    if (j > 0){
        neighbors[Nnb] = ij.previous_y();
        ++Nnb;
    }
    if (j < ij.ny-1){
        neighbors[Nnb] = ij.next_y();
        ++Nnb;
    }

    double angle_ij = angle[ij.flat()];
    double s_ij = sign[ij.flat()];

    xy_t vi(s_ij * std::sin(angle_ij), s_ij * std::cos(angle_ij));

    double cost = -(Nnb * vi.dot(v_global)) / cost_beta;
    for (size_t k=0; k<Nnb; ++k){
        const GridIndex<uint32_t>& nb = neighbors[k];
        s_ij = sign[nb.flat()];
        angle_ij = angle[nb.flat()];
        xy_t vn(s_ij * std::sin(angle_ij), s_ij * std::cos(angle_ij));
        cost -= vi.dot(vn);
    }

    return cost;
}


void unwrap_azimuth_field(
        double* angle, uint32_t nx, uint32_t ny, size_t Nmax,
        double cost_beta
)
{
    /*
     * Tries to unravel the angle and fix discontinuities.
     * 'cost_beta' weights agreement with the global average direction
     * vs. local agreement.
     */


    /* A structure to sort linear indices by associated cost: */
    struct node_cost_t {
        size_t i_lin;
        double cost;

        node_cost_t(size_t i_lin, double cost) : i_lin(i_lin), cost(cost)
        {}

        node_cost_t() = default;

        bool operator<(const node_cost_t& other) const
        {
            return cost < other.cost;
        }
    };

    std::priority_queue<node_cost_t> queue;

    const size_t N = static_cast<size_t>(nx) * ny;
    const size_t age_max = std::max<size_t>(round((0.01*nx)*ny), 4);

    /* Initialize the signs, all positive: */
    std::vector<int8_t> sign(N, 1);

    /* Initial mean velocity: */
    xy_t v_global(global_mean_velocity(angle, sign, nx, ny));

    /* This recomputes the cost and refills the queue. */
    auto initialize_queue_and_compute_cost
    = [nx, ny, cost_beta, angle, &v_global, &sign, &queue]() -> long double
    {
        if (!queue.empty())
            queue = std::priority_queue<node_cost_t>();

        long double new_cost = 0.0;
        for (size_t i=0; i<nx; ++i){
            for (size_t j=0; j<ny; ++j){
                GridIndex<uint32_t> ij(i, j, nx, ny);
                double cost_ij = cost_contribution(
                    ij, angle, sign, v_global, nx, cost_beta
                );
                new_cost += cost_ij;
                sign[ij.flat()] *= -1;
                double cost_delta
                    = cost_contribution(
                        ij, angle, sign, v_global, nx, cost_beta
                    ) - cost_ij;
                sign[ij.flat()] *= -1;
                if (cost_delta < 0)
                    queue.emplace(ij.flat(), cost_delta);
            }
        }

        return new_cost;
    };

    long double cost = initialize_queue_and_compute_cost();

    /*
     * Computing the cost delta for 180° flipping a single index.
     */
    auto proposed_cost_delta =
    [nx, cost_beta, angle, &sign, &v_global]
    (const GridIndex<uint32_t>& id) -> double
    {
        double cost_ij = cost_contribution(
            id, angle, sign, v_global, nx, cost_beta);
        sign[id.flat()] *= -1;
        double cost_delta
            = cost_contribution(
                id, angle, sign, v_global, nx, cost_beta
            ) - cost_ij;
        sign[id.flat()] *= -1;
        return cost_delta;
    };

    size_t l = 0;
    size_t vg_age = 0;
    for (size_t k=0; k<Nmax; ++k){
        /* Increase the age: */
        ++vg_age;

        /* If the global velocity is too old (in terms of sign updates),
         * we wanna recompute it.
         * Also if the queue is empty, recompute and reinitialize the queue
         * so as to possibly squeeze out the last cost changes out: */
        if (vg_age > age_max || queue.empty()){
            v_global = global_mean_velocity(angle, sign, nx, ny);
            vg_age = 0;

            cost = initialize_queue_and_compute_cost();
        }

        /* If now the queue is empty, we cannot further reduce the cost! */
        if (queue.empty())
            break;

        /* Find the element with the largest negative contribution.
         * Ensure that the contribution is actually negative!
         * We might have kept old indices that are already obsolete.
         */
        size_t i_lin = queue.top().i_lin;
        queue.pop();
        GridIndex<uint32_t> current
            = GridIndex<uint32_t>::from_linear(i_lin, nx, ny);
        double cost_delta = proposed_cost_delta(current);


        while (!queue.empty() && cost_delta >= 0.0)
        {
            i_lin = queue.top().i_lin;
            current = GridIndex<uint32_t>::from_linear(i_lin, nx, ny);
            queue.pop();
            cost_delta = proposed_cost_delta(current);
        }

        if (cost_delta >= 0.0)
            /* No success. */
            break;

        /* Accept the flip: */
        sign[current.flat()] *= -1;
        cost += cost_delta;
        l += 1;

        /* Check all the neighbors for potentially adding them */
        std::array<GridIndex<uint32_t>,4> neighbors;
        uint_fast8_t Nnb = 0;
        if (current.i > 0){
            neighbors[0] = current.previous_x();
            ++Nnb;
        }
        if (current.i < nx-1){
            neighbors[Nnb] = current.next_x();
            ++Nnb;
        }
        if (current.j > 0){
            neighbors[Nnb] = current.previous_y();
            ++Nnb;
        }
        if (current.j < ny-1){
            neighbors[Nnb] = current.next_y();
            ++Nnb;
        }

        for (uint_fast8_t k=0; k<Nnb; ++k){
            double cost_delta = proposed_cost_delta(neighbors[k]);
            if (cost_delta < 0){
                queue.emplace(neighbors[k].flat(), cost_delta);
            }
        }
    }

    /* Now perform all the flips: */
    if (l > 0){
        for (size_t k=0; k<N; ++k){
            if (sign[k] == -1){
                double& ak = angle[k];
                if (ak <= 0.0)
                    ak += std::numbers::pi_v<double>;
                else
                    ak -= std::numbers::pi_v<double>;
            }
        }
    }

}

} // end namespace