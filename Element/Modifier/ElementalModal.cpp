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

#include "ElementalModal.h"

ElementalModal::ElementalModal(const unsigned T, const double F, const double D, uvec&& ET)
    : ModifierDynamics(T, std::move(ET))
    , cut_off_freq(F * F)
    , damping(2. * D) {}

int ElementalModal::update_status() {
    suanpan::for_all(element_pool, [&](const weak_ptr<Element>& ele_ptr) {
        const auto t_ptr = ele_ptr.lock();

        if(nullptr == t_ptr || !t_ptr->if_update_viscous() || !t_ptr->allow_modify_viscous()) return;

        if(t_ptr->get_current_mass().empty() || t_ptr->get_current_stiffness().empty()) return;

        cx_vec eigval;
        cx_mat eigvec;
        if(!eig_pair(eigval, eigvec, t_ptr->get_current_stiffness(), t_ptr->get_current_mass())) return;
        cx_mat t_damping(t_ptr->get_total_number(), t_ptr->get_total_number(), fill::zeros);

        const cx_mat theta = t_ptr->get_current_mass() * eigvec;

        for(uword I = 0; I < eigval.n_elem; ++I) if(abs(eigval(I)) < cut_off_freq) t_damping += theta.col(I) * theta.col(I).t() * damping * sqrt(eigval(I)) / dot(theta.col(I), t_ptr->get_current_mass() * theta.col(I));

        access::rw(t_ptr->get_trial_damping_force()) = (access::rw(t_ptr->get_trial_viscous()) = abs(t_damping)) * get_trial_velocity(t_ptr.get());
    });

    return SUANPAN_SUCCESS;
}
