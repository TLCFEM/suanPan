/*******************************************************************************
 * Copyright (C) 2017-2026 Theodore Chang
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

#include "NodeGroup.h"

#include <Domain/DomainBase.h>
#include <Domain/Node.h>

NodeGroup::NodeGroup(const unsigned T, const int D, vec&& R)
    : Group(T)
    , dof(D - 1)
    , rule(std::move(R)) {}

NodeGroup::NodeGroup(const unsigned T, uvec&& R)
    : Group(T, std::move(R))
    , dof(-1) {}

NodeGroup::NodeGroup(const unsigned T, vec&& SN, vec&& EN)
    : Group(T)
    , dof(-2)
    , s_node(std::move(SN))
    , e_node(std::move(EN)) {}

NodeGroup::NodeGroup(const unsigned T, vec&& PL)
    : Group(T)
    , dof(-3)
    , rule(std::move(PL)) {}

void NodeGroup::initialize(const shared_ptr<DomainBase>& D) {
    // generate by direct assignment
    if(-1 == dof) return;

    std::vector<uword> pond;

    if(-2 == dof)
        // generate by line segment between two points
        for(auto& I : D->get_node_pool()) {
            auto& J = I->get_coordinate();

            const auto size = std::min(J.n_elem, std::min(s_node.n_elem, e_node.n_elem));

            if(0 == size) continue;

            const vec s_point(s_node.mem, size), e_point(e_node.mem, size);
            const vec line_b = e_point - s_point;

            const auto dot_bb = dot(line_b, line_b);
            if(std::sqrt(dot_bb) <= datum::eps) continue;

            const vec m_point(J.mem, size);
            const vec line_a = m_point - s_point;

            const auto dot_ab = dot(line_a, line_b);

            if(const auto proj_ab = dot_ab / dot_bb; std::fabs(dot(line_a, line_a) - dot_ab * proj_ab) <= datum::eps * dot_bb && proj_ab >= 0. && proj_ab <= 1.) pond.emplace_back(I->get_tag());
        }
    else if(-3 == dof) {
        // generate by plane fitting
        const auto ref_error = norm(rule);
        for(auto& I : D->get_node_pool()) {
            auto& J = I->get_coordinate();

            const auto size = std::min(J.n_elem, rule.n_elem - 1);

            if(0 == size) continue;

            if(const vec part(rule.mem, size), m_point(J.mem, size); std::fabs(dot(part, m_point) + rule.back()) <= datum::eps * ref_error) pond.emplace_back(I->get_tag());
        }
    }
    else {
        // generate by polynomial curve fitting
        const auto ref_error = norm(rule);
        for(auto& I : D->get_node_pool())
            if(auto& J = I->get_coordinate(); std::cmp_greater(J.n_elem, dof) && std::fabs(as_scalar(polyval(rule, vec{J(dof)}))) <= datum::eps * ref_error) pond.emplace_back(I->get_tag());
    }

    suanpan_sort(pond.begin(), pond.end());

    pool = pond;
}
