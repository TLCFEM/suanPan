/*******************************************************************************
 * Copyright (C) 2017-2022 Theodore Chang
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

double ParticleCollision2D::compute_f(const double distance) const { return distance >= space ? 0. : -alpha * log(distance / space); }

double ParticleCollision2D::compute_df(const double distance) const { return distance >= space ? 0. : -alpha / distance; }

int ParticleCollision2D::process_meta(const shared_ptr<DomainBase>& D, const bool full) {
    auto& W = D->get_factory();

    auto& node_pool = D->get_node_pool();

    list = std::vector<CellList>(node_pool.size());

    decltype(list.size()) counter = 0;
    for(auto& I : node_pool) {
        if(norm(I->get_trial_velocity()) * W->get_incre_time() > space) suanpan_warning("the maximum particle speed seems to be too large, please decrease time increment.\n");
        const auto new_pos = get_position(I);
        list[counter].y = static_cast<int>(floor(new_pos(1) / space));
        list[counter].x = static_cast<int>(floor(new_pos(0) / space));
        list[counter++].tag = I->get_tag();
    }

    suanpan_sort(list.begin(), list.end(), [](const CellList& a, const CellList& b) { return a.x < b.x || a.x == b.x && a.y < b.y; });

    resistance.zeros(W->get_size());

    for(auto I = list.cbegin(); I != list.cend(); ++I)
        for(auto J = I + 1; J != list.cend(); ++J) {
            const auto diff_x = J->x - I->x;
            if(diff_x > 1) break;
            const auto diff_y = J->y - I->y;
            if(diff_x == 1 && diff_y > 1) break;
            if(abs(diff_y) > 1) continue;
            apply_contact(D, D->get<Node>(I->tag), D->get<Node>(J->tag), full);
        }

    return SUANPAN_SUCCESS;
}

ParticleCollision2D::ParticleCollision2D(const unsigned T, const unsigned S, const double G, const double A)
    : ParticleCollision(T, S, 2)
    , space(std::max(G, 1E-4))
    , alpha(A) {}
