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

#include "LineUDL.h"

#include <Domain/Factory.hpp>

LineUDL::LineUDL(const unsigned T, const double L, uvec&& N, std::vector<Node::DOF>&& DT, const unsigned AT, const unsigned D)
    : Load(T, AT, suanpan::translational(D), std::move(DT), L)
    , dimension(D) { target_node = std::move(N); }

int LineUDL::initialize(const shared_ptr<DomainBase>& D) {
    if(SUANPAN_SUCCESS != Load::initialize(D)) return SUANPAN_FAIL;

    if(!validate_node(D)) return SUANPAN_FAIL;

    target_dof = collect_node_dof(D);

    return SUANPAN_SUCCESS;
}

int LineUDL::process(const shared_ptr<DomainBase>& D) {
    if(target_dof.is_empty()) return SUANPAN_SUCCESS;

    trial_load.zeros(D->get_factory()->get_size());

    D->insert_loaded_dof(target_dof);

    mat distribution(dimension, target_node.n_elem, fill::zeros);
    for(uword I{0}, J{1}; J < target_node.n_elem; ++I, ++J) {
        const auto projection = project(abs(D->get<Node>(target_node(J))->initial_position(dimension) - D->get<Node>(target_node(I))->initial_position(dimension)));
        distribution.col(I) += projection;
        distribution.col(J) += projection;
    }

    const auto ref_load = .5 * magnitude * get_amplitude(D);
    if(const auto tag = get_dof_component()[0]; Node::DOF::U1 == tag) trial_load(target_dof) = ref_load * distribution.row(0).t();
    else if(Node::DOF::U2 == tag) trial_load(target_dof) = ref_load * distribution.row(1).t();
    else if(Node::DOF::U3 == tag && 3u == dimension) trial_load(target_dof) = ref_load * distribution.row(2).t();

    return SUANPAN_SUCCESS;
}

LineUDL2D::LineUDL2D(const unsigned T, const double L, uvec&& N, std::vector<Node::DOF>&& DT, const unsigned AT)
    : LineUDL(T, L, std::move(N), std::move(DT), AT, 2u) {}

vec LineUDL2D::project(vec&& increment) const { return reverse(increment); }

LineUDL3D::LineUDL3D(const unsigned T, const double L, uvec&& N, std::vector<Node::DOF>&& DT, const unsigned AT)
    : LineUDL(T, L, std::move(N), std::move(DT), AT, 3u) {}

vec LineUDL3D::project(vec&& increment) const {
    increment.transform([](const double value) { return value * value; });
    return sqrt(vec{increment(1) + increment(2), increment(0) + increment(2), increment(0) + increment(1)});
}
