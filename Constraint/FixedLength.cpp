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

#include "FixedLength.h"

#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Domain/Node.h>

FixedLength::FixedLength(const unsigned T, const unsigned D, uvec&& N)
    : Constraint(T, 0, std::move(N), 2 == D ? uvec{1, 2} : uvec{1, 2, 3}, 1) { set_connected(true); }

int FixedLength::initialize(const shared_ptr<DomainBase>& D) {
    dof_encoding = get_nodal_active_dof(D);

    // need to check if sizes conform since the method does not emit error flag
    if(dof_encoding.n_elem != node_encoding.n_elem * dof_reference.n_elem) return SUANPAN_FAIL;

    coor = resize(D->get<Node>(node_encoding(1))->get_coordinate(), dof_reference.n_elem, 1) - resize(D->get<Node>(node_encoding(0))->get_coordinate(), dof_reference.n_elem, 1);

    set_multiplier_size(0);

    return Constraint::initialize(D);
}

int FixedLength::process(const shared_ptr<DomainBase>& D) {
    auto& W = D->get_factory();

    const auto& n_dof = dof_reference.n_elem;

    const uvec dof_i = dof_encoding.head(n_dof);
    const uvec dof_j = dof_encoding.tail(n_dof);

    const vec t_disp = W->get_trial_displacement()(dof_j) - W->get_trial_displacement()(dof_i);

    if(const auto t_gap = accu(square(coor + t_disp)); min_bound && max_bound) {
        if(0u == num_size && t_gap > min_gap && t_gap < max_gap) return SUANPAN_SUCCESS;

        auxiliary_load = (2. * std::sqrt(t_gap) < std::sqrt(min_gap) + std::sqrt(max_gap) ? min_gap : max_gap) - dot(coor, coor);
    }
    else if(min_bound && !max_bound) {
        if(0u == num_size && t_gap > min_gap) return SUANPAN_SUCCESS;

        auxiliary_load = min_gap - dot(coor, coor);
    }
    else if(!min_bound && max_bound) {
        if(0u == num_size && t_gap < max_gap) return SUANPAN_SUCCESS;

        auxiliary_load = max_gap - dot(coor, coor);
    }

    set_multiplier_size(1);

    auxiliary_stiffness.zeros(W->get_size(), num_size);
    auxiliary_resistance = 0.;
    for(auto I = 0llu; I < n_dof; ++I) {
        auxiliary_stiffness(dof_i(I)) = -(auxiliary_stiffness(dof_j(I)) = 2. * (coor(I) + t_disp(I)));
        auxiliary_resistance += t_disp(I) * (2. * coor(I) + t_disp(I));
    }

    stiffness.zeros(dof_encoding.n_elem, dof_encoding.n_elem);
    const auto t_factor = 2. * trial_lambda(0);
    for(auto I = 0llu; I < n_dof; ++I) stiffness(I + n_dof, I) = stiffness(I, I + n_dof) = -(stiffness(I, I) = stiffness(I + n_dof, I + n_dof) = t_factor);

    resistance = auxiliary_stiffness * trial_lambda;

    return SUANPAN_SUCCESS;
}

void FixedLength::update_status(const vec& i_lambda) { trial_lambda += i_lambda; }

void FixedLength::commit_status() {
    current_lambda = trial_lambda;
    set_multiplier_size(0);
}

void FixedLength::clear_status() {
    current_lambda = trial_lambda.zeros();
    set_multiplier_size(0);
}

void FixedLength::reset_status() {
    trial_lambda = current_lambda;
    set_multiplier_size(0);
}

MinimumGap::MinimumGap(const unsigned T, const unsigned D, const double M, uvec&& N)
    : FixedLength(T, D, std::move(N)) {
    min_bound = true;
    min_gap = M * M;
}

MaximumGap::MaximumGap(const unsigned T, const unsigned D, const double M, uvec&& N)
    : FixedLength(T, D, std::move(N)) {
    max_bound = true;
    max_gap = M * M;
}

Sleeve::Sleeve(const unsigned T, const unsigned D, const double M1, const double M2, uvec&& N)
    : FixedLength(T, D, std::move(N)) {
    min_bound = true;
    max_bound = true;
    min_gap = M1 * M1;
    max_gap = M2 * M2;
}

MaxForce::MaxForce(const unsigned T, const unsigned D, const double MF, uvec&& N)
    : FixedLength(T, D, std::move(N))
    , max_force(MF) {}

int MaxForce::process(const shared_ptr<DomainBase>& D) {
    if(current_flag) {
        // if already exceeded, the constraint is not triggered
        set_multiplier_size(0);
        return SUANPAN_SUCCESS;
    }

    if(SUANPAN_SUCCESS != FixedLength::process(D)) return SUANPAN_FAIL;

    if(0u == num_size) return SUANPAN_SUCCESS;

    vec nodal_resistance(dof_reference.n_elem);
    for(auto I = 0llu; I < nodal_resistance.n_elem; ++I) nodal_resistance(I) = resistance(dof_encoding(I));

    if(norm(nodal_resistance) > max_force) {
        trial_flag = true;
        set_multiplier_size(0);
    }

    return SUANPAN_SUCCESS;
}

void MaxForce::commit_status() {
    current_flag = trial_flag;
    return FixedLength::commit_status();
}

void MaxForce::clear_status() {
    current_flag = trial_flag = false;
    return FixedLength::clear_status();
}

void MaxForce::reset_status() {
    trial_flag = current_flag;
    return FixedLength::reset_status();
}
