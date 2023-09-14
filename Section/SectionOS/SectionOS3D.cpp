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

#include "SectionOS3D.h"
#include <Material/Material.h>

const mat SectionOS3D::weighing_mat = [] {
    mat X(4, 4, fill::zeros);
    X.diag().fill(1. / 15.);
    X(0, 1) = X(1, 0) = X(2, 3) = X(3, 2) = -1. / 60.;
    return X;
}();

SectionOS3D::IntegrationPoint::IntegrationPoint(const double CY, const double CZ, const double CS, const double CN, const double W, unique_ptr<Material>&& M)
    : coor_y(CY)
    , coor_z(CZ)
    , coor_s(CS)
    , coor_n(CN)
    , weight(W)
    , s_material(std::forward<unique_ptr<Material>>(M)) {}

SectionOS3D::SectionOS3D(const unsigned T, const unsigned MT, const double A, vec&& E)
    : Section(T, SectionType::OS3D, MT, A, std::forward<vec>(E)) {}

void SectionOS3D::register_elemental_deformation(const vec& t_deformation) { elemental_deformation = t_deformation; }

/**
 * \brief The deformation is assumed to contain the following.
 *
 *  [0]: u'       \n
 *  [1]: v''      \n
 *  [2]: w''      \n
 *  [3]: f        \n
 *  [4]: f'       \n
 *  [5]: f''      \n
 *  [6]: theta_zi \n
 *  [7]: theta_zj \n
 *  [8]: theta_yi \n
 *  [9]: theta_yj \n
 *
 */
int SectionOS3D::update_trial_status(const vec& t_deformation) {
    if(const vec incre_deformation = (trial_deformation = t_deformation) - current_deformation; norm(incre_deformation) <= datum::eps) return SUANPAN_SUCCESS;

    const auto& up = trial_deformation(0);
    const auto& vpp = trial_deformation(1);
    const auto& wpp = trial_deformation(2);
    const auto& f = trial_deformation(3);
    const auto& fp = trial_deformation(4);
    const auto& fpp = trial_deformation(5);

    const auto rotation = elemental_deformation.subvec(1, 4);
    const auto axial_base = dot(rotation, weighing_mat * rotation);

    trial_stiffness.zeros();
    trial_resistance.zeros();

    for(const auto& I : int_pt) {
        const auto arm_y = I.coor_y - eccentricity(0);
        const auto arm_z = I.coor_z - eccentricity(1);

        const auto axial_strain = axial_base + up - arm_y * vpp - arm_z * wpp + (arm_z * vpp - arm_y * wpp) * f + .5 * (arm_y * arm_y + arm_z * arm_z) * fp * fp + I.coor_s * fpp;
        const auto shear_strain = -2. * I.coor_n * fp;

        if(I.s_material->update_trial_status(vec{axial_strain, shear_strain}) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

        mat de(2, 6, fill::zeros);
        de(0, 0) = 1.;
        de(0, 1) = -arm_y + arm_z * f;
        de(0, 2) = -arm_z - arm_y * f;
        de(0, 3) = arm_z * vpp - arm_y * wpp;
        de(0, 4) = (arm_y * arm_y + arm_z * arm_z) * fp;
        de(0, 5) = I.coor_s;
        de(1, 4) = -2. * I.coor_n;

        trial_resistance += I.weight * de.t() * I.s_material->get_trial_stress();
        trial_stiffness += I.weight * de.t() * I.s_material->get_trial_stiffness() * de;
    }

    return SUANPAN_SUCCESS;
}

int SectionOS3D::clear_status() {
    current_deformation = trial_deformation.zeros();
    current_resistance = trial_resistance.zeros();
    current_stiffness = trial_stiffness = initial_stiffness;
    auto code = 0;
    for(const auto& I : int_pt) code += I.s_material->clear_status();
    return code;
}

int SectionOS3D::commit_status() {
    current_deformation = trial_deformation;
    current_resistance = trial_resistance;
    current_stiffness = trial_stiffness;
    auto code = 0;
    for(const auto& I : int_pt) code += I.s_material->commit_status();
    return code;
}

int SectionOS3D::reset_status() {
    trial_deformation = current_deformation;
    trial_resistance = current_resistance;
    trial_stiffness = current_stiffness;
    auto code = 0;
    for(const auto& I : int_pt) code += I.s_material->reset_status();
    return code;
}
