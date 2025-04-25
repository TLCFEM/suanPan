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

#include "Damper01.h"

#include <Domain/DomainBase.h>
#include <Material/Material1D/Material1D.h>

Damper01::Damper01(const unsigned T, uvec&& NT, const unsigned D, const unsigned DIM)
    : MaterialElement1D(T, d_node, 3 == DIM ? 3 : 2, std::move(NT), uvec{D}, false, 3 == DIM ? std::vector{DOF::U1, DOF::U2, DOF::U3} : std::vector{DOF::U1, DOF::U2})
    , d_dof(3 == DIM ? 3 : 2)
    , IS(3 == d_dof ? uvec{0, 1, 2} : uvec{0, 1})
    , JS(3 == d_dof ? uvec{3, 4, 5} : uvec{2, 3}) { modify_viscous = false; }

int Damper01::initialize(const shared_ptr<DomainBase>& D) {
    damper = D->get<Material>(material_tag(0))->get_copy();

    const mat coord = get_coordinate(d_dof);

    access::rw(direction_cosine) = normalise(coord.row(1) - coord.row(0)).t();

    const auto t_disp = get_current_displacement();
    const auto t_vec = get_current_velocity();

    access::rw(damper->get_current_strain()) = vec{dot(direction_cosine, t_disp(JS) - t_disp(IS))};
    access::rw(damper->get_current_strain_rate()) = vec{dot(direction_cosine, t_vec(JS) - t_vec(IS))};

    initial_viscous.set_size(d_size, d_size);
    initial_viscous(IS, IS) = direction_cosine * damper->get_initial_damping() * direction_cosine.t();
    initial_viscous(IS, JS) = -initial_viscous(IS, IS);
    initial_viscous(JS, JS) = initial_viscous(IS, IS);
    initial_viscous(JS, IS) = initial_viscous(IS, JS);

    trial_viscous = current_viscous = initial_viscous;

    if(!damper->get_initial_stiffness().empty()) {
        initial_stiffness.set_size(d_size, d_size);
        initial_stiffness(IS, IS) = direction_cosine * damper->get_initial_stiffness() * direction_cosine.t();
        initial_stiffness(IS, JS) = -initial_stiffness(IS, IS);
        initial_stiffness(JS, JS) = initial_stiffness(IS, IS);
        initial_stiffness(JS, IS) = initial_stiffness(IS, JS);

        trial_stiffness = current_stiffness = initial_stiffness;
    }

    return SUANPAN_SUCCESS;
}

int Damper01::update_status() {
    if(const auto t_disp = get_trial_displacement(), t_vec = get_trial_velocity(); SUANPAN_SUCCESS != damper->update_trial_status(dot(direction_cosine, t_disp(JS) - t_disp(IS)), dot(direction_cosine, t_vec(JS) - t_vec(IS)))) return SUANPAN_FAIL;

    trial_resistance.set_size(d_size);
    trial_resistance(JS) = direction_cosine * damper->get_trial_stress();
    trial_resistance(IS) = -trial_resistance(JS);

    trial_viscous(IS, IS) = direction_cosine * damper->get_trial_damping() * direction_cosine.t();
    trial_viscous(IS, JS) = -trial_viscous(IS, IS);
    trial_viscous(JS, JS) = trial_viscous(IS, IS);
    trial_viscous(JS, IS) = trial_viscous(IS, JS);

    if(!damper->get_trial_stiffness().empty()) {
        trial_stiffness.set_size(d_size, d_size);
        trial_stiffness(IS, IS) = direction_cosine * damper->get_trial_stiffness() * direction_cosine.t();
        trial_stiffness(IS, JS) = -trial_stiffness(IS, IS);
        trial_stiffness(JS, JS) = trial_stiffness(IS, IS);
        trial_stiffness(JS, IS) = trial_stiffness(IS, JS);
    }

    return SUANPAN_SUCCESS;
}

int Damper01::commit_status() { return damper->commit_status(); }

int Damper01::clear_status() { return damper->clear_status(); }

int Damper01::reset_status() { return damper->reset_status(); }

std::vector<vec> Damper01::record(const OutputType P) { return damper->record(P); }

void Damper01::print() {
    suanpan_info("A viscous damper element using displacement and velocity as basic quantities.\n");
}

int Damper05::update_status() {
    if(const auto t_vec = get_trial_velocity(); SUANPAN_SUCCESS != damper->update_trial_status(0., dot(direction_cosine, t_vec(JS) - t_vec(IS)))) return SUANPAN_FAIL;

    trial_viscous_force.set_size(d_size);
    trial_viscous_force(JS) = direction_cosine * damper->get_trial_stress();
    trial_viscous_force(IS) = -trial_viscous_force(JS);

    trial_viscous(IS, IS) = direction_cosine * damper->get_trial_damping() * direction_cosine.t();
    trial_viscous(IS, JS) = -trial_viscous(IS, IS);
    trial_viscous(JS, JS) = trial_viscous(IS, IS);
    trial_viscous(JS, IS) = trial_viscous(IS, JS);

    return SUANPAN_SUCCESS;
}
