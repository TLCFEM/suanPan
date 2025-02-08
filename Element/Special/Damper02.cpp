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

#include "Damper02.h"
#include <Domain/DomainBase.h>
#include <Material/Material1D/Viscosity/Maxwell.h>

Damper02::Damper02(const unsigned T, uvec&& NT, const unsigned DT, const unsigned ST, const bool UM, const unsigned PC, const double BT, const unsigned DIM)
    : MaterialElement1D(T, d_node, 3 == DIM ? 3 : 2, std::move(NT), {}, false, 3 == DIM ? vector{DOF::U1, DOF::U2, DOF::U3} : vector{DOF::U1, DOF::U2})
    , d_dof(3 == DIM ? 3 : 2)
    , IS(3 == d_dof ? uvec{0, 1, 2} : uvec{0, 1})
    , JS(3 == d_dof ? uvec{3, 4, 5} : uvec{2, 3})
    , device(make_unique<Maxwell>(0, DT, ST, UM, PC, BT)) { modify_viscous = false; }

int Damper02::initialize(const shared_ptr<DomainBase>& D) {
    if(SUANPAN_SUCCESS != device->initialize_base(D) || SUANPAN_SUCCESS != device->initialize(D)) return SUANPAN_FAIL;

    const mat coord = get_coordinate(d_dof);

    access::rw(direction_cosine) = normalise(coord.row(1) - coord.row(0)).t();

    const auto t_disp = get_current_displacement();
    const auto t_vec = get_current_velocity();

    access::rw(device->get_current_strain()) = vec{dot(direction_cosine, t_disp(JS) - t_disp(IS))};
    access::rw(device->get_current_strain_rate()) = vec{dot(direction_cosine, t_vec(JS) - t_vec(IS))};

    initial_viscous.set_size(d_size, d_size);
    initial_viscous(IS, IS) = direction_cosine * device->get_initial_damping() * direction_cosine.t();
    initial_viscous(IS, JS) = -initial_viscous(IS, IS);
    initial_viscous(JS, JS) = initial_viscous(IS, IS);
    initial_viscous(JS, IS) = initial_viscous(IS, JS);

    initial_stiffness.set_size(d_size, d_size);
    initial_stiffness(IS, IS) = direction_cosine * device->get_initial_stiffness() * direction_cosine.t();
    initial_stiffness(IS, JS) = -initial_stiffness(IS, IS);
    initial_stiffness(JS, JS) = initial_stiffness(IS, IS);
    initial_stiffness(JS, IS) = initial_stiffness(IS, JS);

    trial_viscous = current_viscous = initial_viscous;
    trial_stiffness = current_stiffness = initial_stiffness;

    return SUANPAN_SUCCESS;
}

int Damper02::update_status() {
    if(const auto t_disp = get_trial_displacement(), t_vec = get_trial_velocity(); SUANPAN_SUCCESS != device->update_trial_status(dot(direction_cosine, t_disp(JS) - t_disp(IS)), dot(direction_cosine, t_vec(JS) - t_vec(IS)))) return SUANPAN_FAIL;

    trial_resistance.set_size(d_size);
    trial_resistance(JS) = direction_cosine * device->get_trial_stress();
    trial_resistance(IS) = -trial_resistance(JS);

    trial_viscous.set_size(d_size, d_size);
    trial_viscous(IS, IS) = direction_cosine * device->get_trial_damping() * direction_cosine.t();
    trial_viscous(IS, JS) = -trial_viscous(IS, IS);
    trial_viscous(JS, JS) = trial_viscous(IS, IS);
    trial_viscous(JS, IS) = trial_viscous(IS, JS);

    trial_stiffness.set_size(d_size, d_size);
    trial_stiffness(IS, IS) = direction_cosine * device->get_trial_stiffness() * direction_cosine.t();
    trial_stiffness(IS, JS) = -trial_stiffness(IS, IS);
    trial_stiffness(JS, JS) = trial_stiffness(IS, IS);
    trial_stiffness(JS, IS) = trial_stiffness(IS, JS);

    return SUANPAN_SUCCESS;
}

int Damper02::commit_status() { return device->commit_status(); }

int Damper02::clear_status() { return device->clear_status(); }

int Damper02::reset_status() { return device->reset_status(); }

vector<vec> Damper02::record(const OutputType P) { return device->record(P); }

void Damper02::print() {
    suanpan_info("A viscous damper element using displacement and velocity as basic quantities.\n");
}
