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

#include "LeeElementalDamping.h"

LeeElementalDamping::LeeElementalDamping(const unsigned T, const double A, const double B, const double C, const double D, uvec&& ET)
    : Modifier(T, std::forward<uvec>(ET))
    , a(A)
    , b(B)
    , c(C)
    , d(D) {}

int LeeElementalDamping::update_status() {
    suanpan::for_all(element_pool, [&](const weak_ptr<Element>& ele_ptr) {
        if(const auto t_ptr = ele_ptr.lock(); nullptr != t_ptr && t_ptr->if_update_damping()) {
            if(a != 0. && !t_ptr->get_current_mass().empty()) access::rw(t_ptr->get_mass_container()) = a * t_ptr->get_current_mass();

            mat t_stiffness(t_ptr->get_total_number(), t_ptr->get_total_number(), fill::zeros);
            if(b != 0. && !t_ptr->get_current_stiffness().empty()) {
                t_stiffness += b * t_ptr->get_current_stiffness();
                if(t_ptr->is_nlgeom() && !t_ptr->get_current_geometry().empty()) t_stiffness += b * t_ptr->get_current_geometry();
            }
            if(c != 0. && !t_ptr->get_initial_stiffness().empty()) t_stiffness += c * t_ptr->get_initial_stiffness();
            if(d != 0. && !t_ptr->get_trial_stiffness().empty()) {
                t_stiffness += d * t_ptr->get_trial_stiffness();
                if(t_ptr->is_nlgeom() && !t_ptr->get_trial_geometry().empty()) t_stiffness += d * t_ptr->get_trial_geometry();
            }

            access::rw(t_ptr->get_stiffness_container()) = t_stiffness;
        }
    });

    return SUANPAN_SUCCESS;
}
