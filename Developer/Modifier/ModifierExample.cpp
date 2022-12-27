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

#include "ModifierExample.h"
#include <Toolbox/utility.h>

SUANPAN_EXPORT void new_modifierexample(unique_ptr<Modifier>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_modifierexample() requires a valid tag.\n");
        return;
    }

    double a, b;
    if(!get_input(command, a)) {
        suanpan_info("new_modifierexample() needs two valid numbers.\n");
        return;
    }
    if(!get_input(command, b)) {
        suanpan_info("new_modifierexample() needs two valid numbers.\n");
        return;
    }

    std::vector<uword> element_tag;
    unsigned e_tag;
    while(!command.eof()) if(get_input(command, e_tag)) element_tag.emplace_back(e_tag);

    return_obj = make_unique<ModifierExample>(tag, a, b, element_tag);
}

ModifierExample::ModifierExample(const unsigned T, const double A, const double B, uvec&& ET)
    : Modifier(T, std::forward<uvec>(ET))
    , a(A)
    , b(B) {}

int ModifierExample::update_status() {
    suanpan::for_all(element_pool, [&](const weak_ptr<Element>& ele_ptr) {
        if(const auto t_ptr = ele_ptr.lock()) {
            mat t_damping(t_ptr->get_total_number(), t_ptr->get_total_number(), fill::zeros);
            if(a != 0. && !t_ptr->get_current_mass().empty()) t_damping += a * t_ptr->get_current_mass();
            if(b != 0. && !t_ptr->get_current_stiffness().empty()) t_damping += b * t_ptr->get_current_stiffness();
            access::rw(t_ptr->get_trial_damping()) = t_damping;
            access::rw(t_ptr->get_trial_damping_force()) = t_damping * get_trial_velocity(t_ptr.get());
        }
    });

    return SUANPAN_SUCCESS;
}
