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

#include "BilinearElastic1D.h"
#include <Toolbox/utility.h>

BilinearElastic1D::BilinearElastic1D(const unsigned T, const double E, const double Y, const double H, const double R, const double D)
    : DataBilinearElastic1D{fabs(E), fabs(Y), fabs(E) * H, R}
    , Material1D(T, D) {}

int BilinearElastic1D::initialize(const shared_ptr<DomainBase>&) {
    trial_stiffness = current_stiffness = initial_stiffness = elastic_modulus;

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> BilinearElastic1D::get_copy() { return make_unique<BilinearElastic1D>(*this); }

int BilinearElastic1D::update_trial_status(const vec& t_strain) {
    incre_strain = (trial_strain = t_strain) - current_strain;

    if(fabs(incre_strain(0)) <= datum::eps) return SUANPAN_SUCCESS;

    if(const auto abs_strain = fabs(trial_strain(0)); suanpan::approx_equal(radius, 0.)) { abs_strain > yield_strain ? trial_stress = yield_stress + (abs_strain - yield_strain) * (trial_stiffness = hardening_modulus) : trial_stress = abs_strain * (trial_stiffness = elastic_modulus); }
    else {
        const auto normal_strain = std::max(datum::eps, abs_strain / yield_strain);
        const auto factor_a = 1. + pow(normal_strain, radius);
        const auto factor_b = (1. - hardening_modulus / elastic_modulus) * pow(factor_a, -1. / radius);

        trial_stress = (hardening_modulus * yield_strain + factor_b * yield_stress) * normal_strain;
        trial_stiffness = hardening_modulus + elastic_modulus * factor_b / factor_a;
    }

    if(suanpan::sign(trial_strain(0)) < 0.) trial_stress = -trial_stress;

    return SUANPAN_SUCCESS;
}

int BilinearElastic1D::clear_status() {
    current_strain.zeros();
    current_stress.zeros();
    current_stiffness = initial_stiffness;
    return reset_status();
}

int BilinearElastic1D::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_stiffness = trial_stiffness;
    return SUANPAN_SUCCESS;
}

int BilinearElastic1D::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_stiffness = current_stiffness;
    return SUANPAN_SUCCESS;
}

void BilinearElastic1D::print() {
    suanpan_info("A bilinear elastic material model with an elastic modulus of {:.3E} and a hardening ratio of {:.2f}.\n", elastic_modulus, hardening_modulus / elastic_modulus);
    Material1D::print();
}
