/*******************************************************************************
 * Copyright (C) 2017-2025 Theodore Chang
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/

#include "ParticleCollision2D.h"

#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Domain/Node.h>

double ParticleCollision2D::compute_f(const double distance) const { return distance >= space ? 0. : -alpha * log(distance / space); }

double ParticleCollision2D::compute_df(const double distance) const { return distance >= space ? 0. : -alpha / distance; }

int ParticleCollision2D::process_meta(const shared_ptr<DomainBase>& D, const bool full) {
    auto& W = D->get_factory();

    auto& node_pool = D->get_node_pool();

    const auto node_size = node_pool.size();

    list = std::vector<CellList>(node_size);

    suanpan::for_each(node_size, [&](const size_t I) {
        const auto& t_node = node_pool[I];
        if(norm(t_node->get_trial_velocity()) * W->get_incre_time() > space)
            suanpan_warning("The nodal speed seems to be too large.\n");
        const auto new_pos = get_position(t_node);
        list[I].y = static_cast<int>(floor(new_pos(1) / space));
        list[I].x = static_cast<int>(floor(new_pos(0) / space));
        list[I].tag = t_node->get_tag();
    });

    suanpan_sort(list.begin(), list.end(), [](const CellList& a, const CellList& b) { return a.x < b.x || a.x == b.x && a.y < b.y; });

    resistance.zeros(W->get_size());

    suanpan::for_each(node_size, [&](const size_t I) {
        for(auto J = I + 1; J < node_size; ++J) {
            const auto diff_x = list[J].x - list[I].x;
            if(diff_x > 1) break;
            const auto diff_y = list[J].y - list[I].y;
            if(diff_x == 1 && diff_y > 1) break;
            if(std::abs(diff_y) > 1) continue;
            apply_contact(D, D->get<Node>(list[I].tag), D->get<Node>(list[J].tag), full);
        }
    });

    return SUANPAN_SUCCESS;
}

ParticleCollision2D::ParticleCollision2D(const unsigned T, const double G, const double A)
    : ParticleCollision(T, 2)
    , space(std::max(G, 1E-4))
    , alpha(A) {}
