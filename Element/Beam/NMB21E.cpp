/*******************************************************************************
 * Copyright (C) 2017-2024 Theodore Chang
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

#include "NMB21E.h"
#include <Domain/DomainBase.h>
#include <Element/Utility/Orientation.h>
#include <Section/Section.h>

NMB21E::NMB21E(const unsigned T, const unsigned W, uvec&& N, const unsigned S, const bool F)
    : NMB21(T, std::move(N), S, F)
    , a{1 == W ? 1llu : 2llu}
    , b{1 == W ? uvec{0llu, 2llu} : uvec{0llu, 1llu}} {}

int NMB21E::initialize(const shared_ptr<DomainBase>& D) {
    if(SUANPAN_SUCCESS != NMB21::initialize(D)) return SUANPAN_FAIL;

    current_local_deformation = trial_local_deformation = b_trans->to_local_vec(get_current_displacement());

    return SUANPAN_SUCCESS;
}

int NMB21E::update_status() {
    b_trans->update_status();

    vec incre_deformation = -trial_local_deformation(b);
    trial_local_deformation(b) = b_trans->to_local_vec(get_trial_displacement())(b);
    incre_deformation += trial_local_deformation(b);

    const auto& t_stiffness = b_section->get_trial_stiffness();
    const auto& t_resistance = b_section->get_trial_resistance();

    trial_local_deformation(a) -= solve(t_stiffness(a, a), t_resistance(a) * length + t_stiffness(a, b) * incre_deformation);

    if(SUANPAN_SUCCESS != b_section->update_trial_status(trial_local_deformation / length)) return SUANPAN_FAIL;

    mat local_stiffness(3, 3, fill::zeros);
    vec local_resistance(3, fill::zeros);
    local_resistance(b) = t_resistance(b) - t_stiffness(b, a) * solve(t_stiffness(a, a), t_resistance(a));
    local_stiffness(b, b) = t_stiffness(b, b) - t_stiffness(b, a) * solve(t_stiffness(a, a), t_stiffness(a, b));

    trial_resistance = b_trans->to_global_vec(local_resistance);
    trial_stiffness = b_trans->to_global_stiffness_mat(local_stiffness / length);

    if(nlgeom) trial_geometry = b_trans->to_global_geometry_mat(local_resistance);

    return SUANPAN_SUCCESS;
}

int NMB21E::commit_status() {
    current_local_deformation = trial_local_deformation;
    return NMB21::commit_status();
}

int NMB21E::clear_status() {
    current_local_deformation = trial_local_deformation.zeros();
    return NMB21::clear_status();
}

int NMB21E::reset_status() {
    trial_local_deformation = current_local_deformation;
    return NMB21::reset_status();
}
