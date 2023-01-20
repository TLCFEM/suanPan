/*******************************************************************************
 * Copyright (C) 2017-2023 Theodore Chang
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

#include "LinearViscosity.h"

LinearViscosity::LinearViscosity(const unsigned T, const double M, uvec&& ET)
    : Modifier(T, std::forward<uvec>(ET))
    , mu{M} {}

int LinearViscosity::update_status() {
    suanpan::for_all(element_pool, [&](const weak_ptr<Element>& ele_ptr) {
        if(const auto t_ptr = ele_ptr.lock(); nullptr != t_ptr && t_ptr->if_update_damping()) {
            mat t_damping(t_ptr->get_total_number(), t_ptr->get_total_number(), fill::zeros);
            t_damping.diag().fill(mu);

            access::rw(t_ptr->get_trial_damping()) = t_damping;
            access::rw(t_ptr->get_trial_damping_force()) = t_damping * get_trial_velocity(t_ptr.get());
        }
    });

    return SUANPAN_SUCCESS;
}
