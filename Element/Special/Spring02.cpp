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

#include "Spring02.h"
#include <Domain/DomainBase.h>
#include <Material/Material.h>

uvec Spring02::IS{0, 1};
uvec Spring02::JS{2, 3};

Spring02::Spring02(const unsigned T, uvec&& NT, const unsigned MT)
    : MaterialElement1D(T, s_node, s_dof, std::forward<uvec>(NT), uvec{MT}, false, {DOF::U1, DOF::U2}) {}

int Spring02::initialize(const shared_ptr<DomainBase>& D) {
    s_material = suanpan::make_copy(D->get<Material>(material_tag(0)));

    const auto t_coord = get_coordinate(2);

    direction_cosine = (t_coord.row(1) - t_coord.row(0)).t();

    access::rw(length) = norm(direction_cosine);

    if(length <= 1E-10) access::rw(length) = 1.;

    direction_cosine /= length * length;

    initial_stiffness.set_size(s_size, s_size);
    initial_stiffness(IS, IS) = direction_cosine * s_material->get_initial_stiffness() * direction_cosine.t();
    initial_stiffness(IS, JS) = -initial_stiffness(IS, IS);
    initial_stiffness(JS, JS) = initial_stiffness(IS, IS);
    initial_stiffness(JS, IS) = initial_stiffness(IS, JS);
    trial_stiffness = current_stiffness = initial_stiffness;

    return SUANPAN_SUCCESS;
}

int Spring02::update_status() {
    if(const auto t_disp = get_trial_displacement(), t_vec = get_trial_velocity(); SUANPAN_SUCCESS != s_material->update_trial_status(dot(direction_cosine, t_disp(JS) - t_disp(IS)), dot(direction_cosine, t_vec(JS) - t_vec(IS)))) return SUANPAN_FAIL;

    trial_stiffness.set_size(s_size, s_size);
    trial_stiffness(IS, IS) = direction_cosine * s_material->get_trial_stiffness() * direction_cosine.t();
    trial_stiffness(IS, JS) = -trial_stiffness(IS, IS);
    trial_stiffness(JS, JS) = trial_stiffness(IS, IS);
    trial_stiffness(JS, IS) = trial_stiffness(IS, JS);

    trial_resistance.set_size(s_size);
    trial_resistance(JS) = direction_cosine * s_material->get_trial_stress();
    trial_resistance(IS) = -trial_resistance(JS);

    return SUANPAN_SUCCESS;
}

int Spring02::commit_status() { return s_material->commit_status(); }

int Spring02::clear_status() { return s_material->clear_status(); }

int Spring02::reset_status() { return s_material->reset_status(); }

vector<vec> Spring02::record(const OutputType P) { return s_material->record(P); }

void Spring02::print() {
    suanpan_info("A spring element that uses strain as basic quantity. The material model used shall be based on displacement--force relationship.\n");
    if(!is_initialized()) return;
    suanpan_info("Material:\n");
    s_material->print();
}
