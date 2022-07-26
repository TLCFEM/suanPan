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

#include "NMB21E.h"
#include <Domain/DomainBase.h>
#include <Element/Utility/Orientation.h>
#include <Section/Section.h>

constexpr unsigned NMB21E::max_iteration = 20;
constexpr double NMB21E::tolerance = 1E-13;

NMB21E::NMB21E(const unsigned T, const unsigned W, uvec&& N, const unsigned S, const bool F)
    : NMB21(T, std::forward<uvec>(N), S, F)
    , which(1 == W ? 1 : 2) {}

int NMB21E::update_status() {
    b_trans->update_status();

    auto local_deformation = b_trans->to_local_vec(get_trial_displacement());
    local_deformation(a) += trial_rotation;

    auto counter = 0u;
    while(true) {
        if(++counter > max_iteration) {
            suanpan_error("NMB21E element %u fails to converge.\n", get_tag());
            return SUANPAN_FAIL;
        }

        if(SUANPAN_SUCCESS != b_section->update_trial_status(local_deformation / length)) return SUANPAN_FAIL;

        mat local_stiffness = b_section->get_trial_stiffness() / length;
        vec local_resistance = b_section->get_trial_resistance();

        const auto error = norm(local_resistance(a));
        suanpan_extra_debug("NMB21E local iteration error: %.4E.\n", error);

        if(error < tolerance) {
            const mat t_mat = local_stiffness(b, b) - local_stiffness(b, a) * solve(local_stiffness(a, a), local_stiffness(a, b));

            local_stiffness.zeros();
            local_stiffness(b, b) = t_mat;

            trial_resistance = b_trans->to_global_vec(local_resistance);
            trial_stiffness = b_trans->to_global_stiffness_mat(local_stiffness);

            if(nlgeom) trial_geometry = b_trans->to_global_geometry_mat(local_resistance);

            return SUANPAN_SUCCESS;
        }

        const vec incre = solve(local_stiffness(a, a), local_resistance(a));
        local_deformation(a) -= incre;
        trial_rotation -= incre;
    }
}

int NMB21E::commit_status() {
    current_rotation = trial_rotation;
    return NMB21::commit_status();
}

int NMB21E::clear_status() {
    current_rotation = trial_rotation.zeros();
    return NMB21::clear_status();
}

int NMB21E::reset_status() {
    trial_rotation = current_rotation;
    return NMB21::reset_status();
}
