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
    : Element(T, D, ET, NT, translational(D))
    , dimension(D)
    , host_tag(ET)
    , alpha(P) {}

int Embedded::initialize(const shared_ptr<DomainBase>& D) {
    host_element = D->get<Element>(host_tag);

    vec normalised_coor(dimension, fill::zeros);

    if(!host_element->is_active() || host_element->compute_shape_function(normalised_coor, 0).empty()) return SUANPAN_FAIL;

    access::rw(host_size) = host_element->get_node_number();

    idx.clear();
    idx.reserve(dimension);
    idx.emplace_back(linspace<uvec>(dimension, dimension * static_cast<uword>(host_size), host_size));
    for(auto I = 1u; I < dimension; ++I) idx.emplace_back(idx.back() + 1);

    const auto temp_coor = get_coordinate(dimension);
    const mat element_coor = temp_coor.tail_rows(host_size);
    const vec node_coor = temp_coor.row(0).t();

    auto counter = 0u;
    while(++counter <= max_iteration) {
        const vec incre = solve((host_element->compute_shape_function(normalised_coor, 1) * element_coor).t(), node_coor - ((shape = host_element->compute_shape_function(normalised_coor, 0)) * element_coor).t());
        if(suanpan::inf_norm(incre) < 1E-14) break;
        normalised_coor += incre;
    }

    return max_iteration == counter ? SUANPAN_FAIL : SUANPAN_SUCCESS;
}

int Embedded::update_status() {
    const mat t_disp = reshape(get_trial_displacement(), dimension, get_node_number());
    const vec reaction = alpha * (t_disp.col(0) - t_disp.tail_cols(host_size) * shape.t());

    trial_resistance.zeros(get_total_number());
    trial_stiffness.zeros(get_total_number(), get_total_number());

    trial_resistance.head(dimension) = -reaction;

    mat t_shape = alpha * shape.t();

    for(auto I = 0u; I < dimension; ++I) {
        trial_resistance(idx[I]) = shape.t() * reaction(I);
        trial_stiffness(I, I) = -alpha;
        trial_stiffness(uvec{I}, idx[I]) = t_shape.t();
        trial_stiffness(idx[I], uvec{I}) = t_shape;
    }

    t_shape *= -shape;

    for(auto I = 0u; I < dimension; ++I) trial_stiffness(idx[I], idx[I]) = t_shape;

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
