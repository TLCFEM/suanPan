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

#include "Bilinear1D.h"

Bilinear1D::Bilinear1D(const unsigned T, const double E, const double Y, const double H, const double B, const double R)
    : DataBilinear1D{fabs(E), fabs(Y), fabs(B), fabs(B * E) * H / (1. - H), fabs((1. - B) * E) * H / (1. - H)}
    , Material1D(T, R) {}

int Bilinear1D::initialize(const shared_ptr<DomainBase>&) {
    trial_stiffness = current_stiffness = initial_stiffness = elastic_modulus;

    initialize_history(2);

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> Bilinear1D::get_copy() { return make_unique<Bilinear1D>(*this); }

int Bilinear1D::update_trial_status(const vec& t_strain) {
    incre_strain = (trial_strain = t_strain) - current_strain;

    if(fabs(incre_strain(0)) <= datum::eps) return SUANPAN_SUCCESS;

    trial_stress = current_stress + (trial_stiffness = initial_stiffness) * incre_strain;

    trial_history = current_history;
    auto& back_stress = trial_history(0);
    auto& plastic_strain = trial_history(1);

    const auto shifted_stress = trial_stress(0) - back_stress;

    const auto yield_surf = yield_stress + isotropic_modulus * plastic_strain;

    if(const auto yield_func = fabs(shifted_stress) - std::max(0., yield_surf); yield_func >= 0.) {
        const auto dkdh = kinematic_modulus + (yield_surf > 0. ? isotropic_modulus : 0.);
        auto incre_plastic_strain = yield_func / (elastic_modulus + dkdh);
        plastic_strain += incre_plastic_strain;
        if(shifted_stress < 0.) incre_plastic_strain = -incre_plastic_strain;
        back_stress += kinematic_modulus * incre_plastic_strain;
        trial_stress -= elastic_modulus * incre_plastic_strain;
        trial_stiffness *= dkdh / (elastic_modulus + dkdh);
    }

    return SUANPAN_SUCCESS;
}

int Bilinear1D::clear_status() {
    current_strain.zeros();
    current_stress.zeros();
    current_history = initial_history;
    current_stiffness = initial_stiffness;
    return reset_status();
}

int Bilinear1D::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_history = trial_history;
    current_stiffness = trial_stiffness;
    return SUANPAN_SUCCESS;
}

int Bilinear1D::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_history = current_history;
    trial_stiffness = current_stiffness;
    return SUANPAN_SUCCESS;
}

void Bilinear1D::print() {
    suanpan_info("A uniaxial bilinear hardening material using J2 plasticity and associated flow rule.\n");
    Material1D::print();
}
