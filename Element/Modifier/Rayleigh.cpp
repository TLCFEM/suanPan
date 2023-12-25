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

#include "Rayleigh.h"

Rayleigh::Rayleigh(const unsigned T, const double A, const double B, const double C, const double D, uvec&& ET)
    : ModifierDynamics(T, std::move(ET))
    , a(A)
    , b(B)
    , c(C)
    , d(D) {}

int Rayleigh::update_status() {
    suanpan::for_all(element_pool, [&](const weak_ptr<Element>& ele_ptr) {
        const auto t_ptr = ele_ptr.lock();

        if(nullptr == t_ptr || !t_ptr->if_update_viscous() || !t_ptr->allow_modify_viscous()) return;

        mat t_damping(t_ptr->get_total_number(), t_ptr->get_total_number(), fill::zeros);
        if(a != 0. && !t_ptr->get_current_mass().empty()) t_damping += a * t_ptr->get_current_mass();
        if(b != 0. && !t_ptr->get_current_stiffness().empty()) {
            t_damping += b * t_ptr->get_current_stiffness();
            if(t_ptr->is_nlgeom() && !t_ptr->get_current_geometry().empty()) t_damping += b * t_ptr->get_current_geometry();
        }
        if(c != 0. && !t_ptr->get_initial_stiffness().empty()) t_damping += c * t_ptr->get_initial_stiffness();
        if(d != 0. && !t_ptr->get_trial_stiffness().empty()) {
            t_damping += d * t_ptr->get_trial_stiffness();
            if(t_ptr->is_nlgeom() && !t_ptr->get_trial_geometry().empty()) t_damping += d * t_ptr->get_trial_geometry();
        }

        access::rw(t_ptr->get_trial_viscous()) = t_damping;

        access::rw(t_ptr->get_trial_damping_force()) = t_damping * get_trial_velocity(t_ptr.get());
    });

    return SUANPAN_SUCCESS;
}
