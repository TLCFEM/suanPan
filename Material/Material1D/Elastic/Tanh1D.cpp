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

#include "Tanh1D.h"

Tanh1D::Tanh1D(const unsigned T, const double E, const double R)
    : DataTanh1D{E}
    , Material1D(T, R) {}

int Tanh1D::initialize(const shared_ptr<DomainBase>&) {
    trial_stiffness = current_stiffness = initial_stiffness = elastic_modulus;

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> Tanh1D::get_copy() { return std::make_unique<Tanh1D>(*this); }

int Tanh1D::update_trial_status(const vec& t_strain) {
    trial_stress = elastic_modulus * arma::tanh(trial_strain = t_strain);
    trial_stiffness = elastic_modulus * arma::pow(arma::cosh(trial_strain), -2.);
    return SUANPAN_SUCCESS;
}

int Tanh1D::clear_status() {
    current_strain = trial_strain.zeros();
    current_stress = trial_stress.zeros();
    current_stiffness = trial_stiffness = elastic_modulus;
    return 0;
}

int Tanh1D::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_stiffness = trial_stiffness;
    return 0;
}

int Tanh1D::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_stiffness = current_stiffness;
    return 0;
}

void Tanh1D::print() {
    suanpan_info("A uniaxial nonlinear elastic material using tanh function with an elastic modulus of {:.4E}.\n", elastic_modulus);
    Material1D::print();
}
