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

#include "NodeLine.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Domain/NodeHelper.hpp>

const mat NodeLine::rotation{{0., -1.}, {1., 0.}};

std::vector<vec> NodeLine::get_position(const shared_ptr<DomainBase>& D) {
    std::vector<vec> position;
    position.reserve(node_encoding.n_elem);

    for(const auto I : node_encoding) position.emplace_back(get_trial_position<DOF::U1, DOF::U2>(D->get<Node>(I)));

    return position;
}

NodeLine::NodeLine(const unsigned T, const unsigned S, const unsigned A, uvec&& N)
    : Constraint(T, S, A, std::move(N), uvec{1, 2}, 1) { set_connected(true); }

int NodeLine::initialize(const shared_ptr<DomainBase>& D) {
    dof_encoding = get_nodal_active_dof(D);

    // need to check if sizes conform since the method does not emit error flag
    if(dof_encoding.n_elem != node_encoding.n_elem * dof_reference.n_elem) return SUANPAN_FAIL;

    set_multiplier_size(0);

    return Constraint::initialize(D);
}

int NodeLine::process(const shared_ptr<DomainBase>& D) {
    const auto node = get_position(D);

    const vec axis = node[1] - node[0];
    const vec outer_normal = rotation * axis;
    const vec position = node[2] - node[0];

    const auto pen = dot(position, outer_normal);

    if(const auto t = dot(position, axis); 0 == num_size && (pen > 0. || t < 0. || t > norm(axis))) return SUANPAN_SUCCESS;

    set_multiplier_size(1);

    const span span_i(0, 1), span_j(2, 3), span_k(4, 5);

    auxiliary_stiffness.zeros(D->get_factory()->get_size(), num_size);
    auxiliary_resistance = pen;

    const rowvec dpdi = -position.t() * rotation - outer_normal.t();
    const rowvec dpdj = position.t() * rotation;

    for(auto I = 0, J = 2, K = 4; I < 2; ++I, ++J, ++K) {
        auxiliary_stiffness(dof_encoding(I)) = dpdi(I);
        auxiliary_stiffness(dof_encoding(J)) = dpdj(I);
        auxiliary_stiffness(dof_encoding(K)) = outer_normal(I);
    }

    const mat factor = trial_lambda(0) * rotation;

    stiffness.zeros(dof_encoding.n_elem, dof_encoding.n_elem);
    stiffness(span_i, span_i) = factor.t() + factor;
    stiffness(span_i, span_j) = stiffness(span_k, span_i) = -(stiffness(span_k, span_j) = factor);
    stiffness(span_i, span_k) = stiffness(span_j, span_i) = -(stiffness(span_j, span_k) = factor.t());

    resistance = auxiliary_stiffness * trial_lambda;

    return SUANPAN_SUCCESS;
}

void NodeLine::update_status(const vec& i_lambda) { trial_lambda += i_lambda; }

void NodeLine::commit_status() {
    current_lambda = trial_lambda;
    set_multiplier_size(0);
}

void NodeLine::clear_status() {
    current_lambda = trial_lambda.zeros();
    set_multiplier_size(0);
}

void NodeLine::reset_status() {
    trial_lambda = current_lambda;
    set_multiplier_size(0);
}
