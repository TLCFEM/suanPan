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

    // generate by two points
    if(-2 == dof)
        for(auto& I : D->get_node_pool()) {
            auto& J = I->get_coordinate();

            const auto size = std::min(J.n_elem, s_node.n_elem);

            if(0 == size) continue;

            const vec s_point(s_node.mem, size);
            const vec e_point(e_node.mem, size);
            const vec m_point(J.mem, size);

            if(const auto denom = e_point(0) - s_point(0); std::fabs(denom) >= 1E-8) {
                if(const auto T = (m_point(0) - s_point(0)) / denom; T >= 0. && T <= 1.) {
                    // on inclined line
                    auto flag = true;
                    for(uword K = 1; K < size; ++K)
                        if(std::fabs((T * ((e_point(K) - s_point(K))) - m_point(K) + s_point(K))) > 1E-8) {
                            flag = false;
                            break;
                        }
                    if(flag) pond.emplace_back(I->get_tag());
                }
            }
            else if(std::fabs(m_point(0) - .5 * (e_point(0) + s_point(0))) < 1E-8) {
                // on vertical line
                auto flag = true;
                for(uword K = 1; K < size; ++K)
                    if(const auto T = (m_point(K) - s_point(K)) / (e_point(K) - s_point(K)); T < 0. || T > 1.) {
                        flag = false;
                        break;
                    }
                if(flag) pond.emplace_back(I->get_tag());
            }
        }
    else if(-3 == dof)
        for(auto& I : D->get_node_pool()) {
            auto& J = I->get_coordinate();

            const auto size = std::min(J.n_elem, rule.n_elem - 1);

            if(0 == size) continue;

            if(const vec part(rule.mem, size), m_point(J.mem, size); std::fabs(dot(part, m_point) + rule.back()) <= 1E-8) pond.emplace_back(I->get_tag());
        }
    else
        // generate by polynomial curve fitting
        for(auto& I : D->get_node_pool())
            if(auto& J = I->get_coordinate(); static_cast<int>(J.n_elem) > dof && std::fabs(as_scalar(polyval(rule, vec{J(dof)}))) <= 1E-12) pond.emplace_back(I->get_tag());

    suanpan_sort(pond.begin(), pond.end());

    pool = pond;
}
