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

#include "B21E.h"
#include <Domain/DomainBase.h>
#include <Element/Utility/Orientation.h>
#include <Section/Section.h>

B21E::B21E(const unsigned T, const unsigned W, uvec&& N, const unsigned S, const unsigned P, const bool F)
    : B21(T, std::move(N), S, P, F)
    , a{1u == W ? 1llu : 2llu}
    , b{1u == W ? uvec{0llu, 2llu} : uvec{0llu, 1llu}} {}

int B21E::update_status() {
    b_trans->update_status();

    auto local_deformation = b_trans->to_local_vec(get_trial_displacement());
    local_deformation(a) += trial_rotation;

    mat local_stiffness(3, 3);
    vec local_resistance(3);

    auto counter = 0u;
    auto ref_error = 1.;
    while(true) {
        if(max_iteration == ++counter) {
            suanpan_error("Element {} fails to converge.\n", get_tag());
            return SUANPAN_FAIL;
        }

        local_stiffness.zeros();
        local_resistance.zeros();

        for(const auto& I : int_pt) {
            if(SUANPAN_SUCCESS != I.b_section->update_trial_status(I.strain_mat * local_deformation / length)) return SUANPAN_FAIL;
            local_stiffness += I.strain_mat.t() * I.b_section->get_trial_stiffness() * I.strain_mat * I.weight / length;
            local_resistance += I.strain_mat.t() * I.b_section->get_trial_resistance() * I.weight;
        }

        const vec incre = solve(local_stiffness(a, a), local_resistance(a));
        const auto error = inf_norm(incre);
        if(1u == counter) ref_error = error;
        suanpan_debug("Local iteration error: {:.5E}.\n", error);

        if(error < tolerance * ref_error || inf_norm(local_resistance(a)) < tolerance) {
            const mat t_mat = local_stiffness(b, b) - local_stiffness(b, a) * solve(local_stiffness(a, a), local_stiffness(a, b));

            local_stiffness.zeros();
            local_stiffness(b, b) = t_mat;

            trial_resistance = b_trans->to_global_vec(local_resistance);
            trial_stiffness = b_trans->to_global_stiffness_mat(local_stiffness);

            if(nlgeom) trial_geometry = b_trans->to_global_geometry_mat(local_resistance);

            return SUANPAN_SUCCESS;
        }

        local_deformation(a) -= incre;
        trial_rotation -= incre;
    }
}

int B21E::commit_status() {
    current_rotation = trial_rotation;
    return B21::commit_status();
}

int B21E::clear_status() {
    current_rotation = trial_rotation.zeros();
    return B21::clear_status();
}

int B21E::reset_status() {
    trial_rotation = current_rotation;
    return B21::reset_status();
}
