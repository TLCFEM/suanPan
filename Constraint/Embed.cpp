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

#include "Embed.h"

#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Domain/Node.h>
#include <Element/Element.h>

Embed::Embed(const unsigned T, const unsigned ET, const unsigned NT, const unsigned D)
    : Constraint(T, 0, {NT}, {}, D)
    , element_tag(ET) {}

int Embed::initialize(const shared_ptr<DomainBase>& D) {
    auto& t_node = D->get<Node>(node_encoding(0));
    auto& t_element = D->get<Element>(element_tag);

    if(nullptr == t_node || nullptr == t_element || !t_node->is_active() || !t_element->is_active() || t_element->compute_shape_function(zeros(num_size, 1), 0).is_empty()) return SUANPAN_FAIL;

    const auto t_coor = get_coordinate(t_element.get(), num_size);

    const vec n_coor = t_node->get_coordinate().head(num_size);

    vec t_para = zeros(num_size);

    rowvec n;

    auto counter = 0u;
    while(true) {
        if(max_iteration == ++counter) return SUANPAN_FAIL;

        const vec incre = solve((t_element->compute_shape_function(t_para, 1) * t_coor).t(), n_coor - ((n = t_element->compute_shape_function(t_para, 0)) * t_coor).t());
        if(inf_norm(incre) < 1E-14) break;
        t_para += incre;
    }

    auto& n_dof = t_node->get_reordered_dof();
    auto& e_dof = t_element->get_dof_encoding();

    auxiliary_stiffness.zeros(D->get_factory()->get_size(), num_size);

    for(auto K = 0u; K < num_size; ++K) {
        auxiliary_stiffness(n_dof(K), K) = -1.;
        for(uword I = 0, J = K; I < n.n_elem; ++I, J += num_size) auxiliary_stiffness(e_dof(J), K) = n(I);
    }

    return Constraint::initialize(D);
}

int Embed::process(const shared_ptr<DomainBase>& D) {
    auxiliary_resistance = auxiliary_stiffness.t() * D->get_factory()->get_trial_displacement();

    return SUANPAN_SUCCESS;
}

Embed2D::Embed2D(const unsigned T, const unsigned ET, const unsigned NT)
    : Embed(T, ET, NT, 2) {}

Embed3D::Embed3D(const unsigned T, const unsigned ET, const unsigned NT)
    : Embed(T, ET, NT, 3) {}
