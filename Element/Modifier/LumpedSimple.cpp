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

#include "LumpedSimple.h"
#include <Element/Utility/MatrixModifier.hpp>

int LumpedSimple::update_status() {
    suanpan::for_all(element_pool, [&](const weak_ptr<Element>& ele_ptr) {
        const auto t_ptr = ele_ptr.lock();

        if(nullptr == t_ptr || !t_ptr->if_update_mass() || !t_ptr->allow_modify_mass()) return;

        if(!t_ptr->get_trial_mass().empty()) suanpan::mass::lumped_simple::apply(access::rw(t_ptr->get_trial_mass()));
    });

    return SUANPAN_SUCCESS;
}
