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

#include "NodeFacet.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Domain/Node.h>
#include <Toolbox/tensorToolbox.h>

std::vector<vec> NodeFacet::get_position(const shared_ptr<DomainBase>& D) {
    std::vector position(node_encoding.n_elem, vec());

    for(auto I = 0llu; I < node_encoding.n_elem; ++I) {
        auto& t_node = D->get<Node>(node_encoding(I));
        position[I] = resize(t_node->get_coordinate(), dof_reference.n_elem, 1) + resize(t_node->get_trial_displacement(), dof_reference.n_elem, 1);
    }

    return position;
}

NodeFacet::NodeFacet(const unsigned T, const unsigned S, const unsigned A, uvec&& N)
    : Constraint(T, S, A, std::forward<uvec>(N), uvec{1, 2, 3}, 1) { set_connected(true); }

int NodeFacet::initialize(const shared_ptr<DomainBase>& D) {
    dof_encoding = get_nodal_active_dof(D);

    if(dof_encoding.n_elem != node_encoding.n_elem * dof_reference.n_elem) return SUANPAN_FAIL;

    set_multiplier_size(0);

    return Constraint::initialize(D);
}

int NodeFacet::process(const shared_ptr<DomainBase>& D) {
    const auto node = get_position(D);

    constexpr auto i = 0, j = 1, k = 2, s = 3;

    const vec edge_i = node[k] - node[j];
    const vec edge_j = node[i] - node[k];
    const vec edge_k = node[j] - node[i];
    const vec outer_normal = cross(edge_j, edge_k);
    vec s_i = node[s] - node[i];

    const auto pen = dot(s_i, outer_normal);

    if(0 == num_size && (pen > 0. || dot(node[s] - node[j], cross(edge_i, outer_normal)) > 0. || dot(node[s] - node[k], cross(edge_j, outer_normal)) > 0. || dot(s_i, cross(edge_k, outer_normal)) > 0.)) return SUANPAN_SUCCESS;

    set_multiplier_size(1);

    const span span_i(0, 2), span_j(3, 5), span_k(6, 8), span_s(9, 11);

    auxiliary_stiffness.zeros(D->get_factory()->get_size(), num_size);
    auxiliary_resistance = pen;

    auto skew_edge_i = transform::skew_symm(edge_i);
    auto skew_edge_j = transform::skew_symm(edge_j);
    auto skew_edge_k = transform::skew_symm(edge_k);

    const rowvec dpdi = s_i.t() * skew_edge_i - outer_normal.t();
    const rowvec dpdj = s_i.t() * skew_edge_j;
    const rowvec dpdk = s_i.t() * skew_edge_k;

    for(auto I = 0, J = 3, K = 6, L = 9; I < 3; ++I, ++J, ++K, ++L) {
        auxiliary_stiffness(dof_encoding(I)) = dpdi(I);
        auxiliary_stiffness(dof_encoding(J)) = dpdj(I);
        auxiliary_stiffness(dof_encoding(K)) = dpdk(I);
        auxiliary_stiffness(dof_encoding(L)) = outer_normal(I);
    }

    skew_edge_i *= trial_lambda(0);
    skew_edge_j *= trial_lambda(0);
    skew_edge_k *= trial_lambda(0);

    const auto skew_s_i = transform::skew_symm(s_i *= trial_lambda(0));

    stiffness.zeros(dof_encoding.n_elem, dof_encoding.n_elem);
    stiffness(span_i, span_j) = -(stiffness(span_j, span_i) = skew_s_i + skew_edge_j);
    stiffness(span_i, span_k) = -(stiffness(span_k, span_i) = skew_edge_k - skew_s_i);
    stiffness(span_j, span_k) = -(stiffness(span_k, span_j) = skew_s_i);
    stiffness(span_i, span_s) = -(stiffness(span_s, span_i) = skew_edge_i);
    stiffness(span_j, span_s) = -(stiffness(span_s, span_j) = skew_edge_j);
    stiffness(span_k, span_s) = -(stiffness(span_s, span_k) = skew_edge_k);

    resistance = auxiliary_stiffness * trial_lambda;

    return SUANPAN_SUCCESS;
}

void NodeFacet::update_status(const vec& i_lambda) { trial_lambda += i_lambda; }

void NodeFacet::commit_status() {
    current_lambda = trial_lambda;
    set_multiplier_size(0);
}

void NodeFacet::clear_status() {
    current_lambda = trial_lambda.zeros();
    set_multiplier_size(0);
}

void NodeFacet::reset_status() {
    trial_lambda = current_lambda;
    set_multiplier_size(0);
}
