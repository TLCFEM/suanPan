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

#include "Embedded.h"

#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>

Embedded::Embedded(const unsigned T, const unsigned ET, const unsigned NT, const unsigned D, const double P)
    : Element(T, D, ET, NT)
    , e_dof(D)
    , host_tag(ET)
    , alpha(P) {}

int Embedded::initialize(const shared_ptr<DomainBase>& D) {
    host_element = D->get<Element>(host_tag);

    if(nullptr == host_element || !host_element->is_active() || host_element->compute_shape_function(zeros(e_dof, 1), 0).empty()) return SUANPAN_FAIL;

    access::rw(host_size) = host_element->get_node_number();

    idx.clear();
    idx.reserve(e_dof);
    idx.emplace_back(linspace<uvec>(e_dof, e_dof * static_cast<uword>(host_size), host_size));
    for(auto I = 1u; I < e_dof; ++I) idx.emplace_back(idx.back() + 1);

    const auto o_coor = get_coordinate(e_dof);
    const mat t_coor = o_coor.tail_rows(host_size);
    const vec n_coor = o_coor.row(0).t();

    vec t_para = zeros(e_dof);

    auto& n = access::rw(iso_n);

    auto counter = 0u;
    while(++counter <= max_iteration) {
        const vec incre = solve((host_element->compute_shape_function(t_para, 1) * t_coor).t(), n_coor - ((n = host_element->compute_shape_function(t_para, 0)) * t_coor).t());
        if(inf_norm(incre) < 1E-14) break;
        t_para += incre;
    }

    return max_iteration == counter ? SUANPAN_FAIL : SUANPAN_SUCCESS;
}

int Embedded::update_status() {
    const auto t_disp = get_trial_displacement();

    const vec reaction = alpha * (t_disp.head(e_dof) - reshape(t_disp.tail(t_disp.n_elem - e_dof), e_dof, host_size) * iso_n.t());

    trial_resistance.zeros(get_total_number());
    trial_stiffness.zeros(get_total_number(), get_total_number());

    trial_resistance.head(e_dof) = -reaction;

    mat t_n = alpha * iso_n.t();

    for(auto I = 0u; I < e_dof; ++I) {
        trial_resistance(idx[I]) = iso_n.t() * reaction(I);
        trial_stiffness(I, I) = -alpha;
        trial_stiffness(uvec{I}, idx[I]) = t_n.t();
        trial_stiffness(idx[I], uvec{I}) = t_n;
    }

    t_n *= -iso_n;

    for(auto I = 0u; I < e_dof; ++I) trial_stiffness(idx[I], idx[I]) = t_n;

    return SUANPAN_SUCCESS;
}

int Embedded::clear_status() {
    trial_resistance.reset();
    trial_stiffness.reset();

    current_resistance = trial_resistance;
    current_stiffness = trial_stiffness;

    return SUANPAN_SUCCESS;
}

int Embedded::commit_status() {
    current_resistance = trial_resistance;
    current_stiffness = trial_stiffness;

    return SUANPAN_SUCCESS;
}

int Embedded::reset_status() {
    trial_resistance = current_resistance;
    trial_stiffness = current_stiffness;

    return SUANPAN_SUCCESS;
}

Embedded2D::Embedded2D(const unsigned T, const unsigned ET, const unsigned NT, const double P)
    : Embedded(T, ET, NT, 2, P) {}

Embedded3D::Embedded3D(const unsigned T, const unsigned ET, const unsigned NT, const double P)
    : Embedded(T, ET, NT, 3, P) {}
