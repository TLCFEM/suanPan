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

#include "LumpedScale.h"

int LumpedScale::update_status() {
    suanpan::for_all(element_pool, [&](const weak_ptr<Element>& ele_ptr) {
        const auto t_ptr = ele_ptr.lock();

        if(nullptr == t_ptr || !t_ptr->if_update_mass() || !t_ptr->allow_modify_mass()) return;

        const auto& t_mass = t_ptr->get_trial_mass();

        if(t_mass.empty()) return;

        const auto num_dof = t_ptr->get_dof_number();
        const auto num_node = t_ptr->get_node_number();
        const auto num_size = num_dof * num_node;
        vec n_mass(num_size, fill::zeros);
        for(unsigned I = 0; I < num_dof; ++I) {
            auto mass_a = 0., mass_b = 0.;
            for(auto J = I; J < num_size; J += num_dof) {
                mass_a += sum(t_mass.col(J));
                mass_b += t_mass(J, J);
            }
            if(mass_b != 0.) {
                const auto factor = mass_a / mass_b;
                for(auto J = I; J < num_size; J += num_dof) n_mass(J) = t_mass(J, J) * factor;
            }
        }

        access::rw(t_ptr->get_trial_mass()) = diagmat(n_mass);
    });

    return SUANPAN_SUCCESS;
}
