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
#include <Recorder/OutputType.h>

const mat SectionOS3D::weighing_mat = [] {
    mat X(4, 4, fill::zeros);
    X.diag().fill(1. / 15.);
    X(0, 1) = X(1, 0) = X(2, 3) = X(3, 2) = -1. / 60.;
    return X;
}();

SectionOS3D::IntegrationPoint::IntegrationPoint(const double CY, const double CZ, const double CS, const double PY, const double PZ, const double W, unique_ptr<Material>&& M)
    : coor_y(CY)
    , coor_z(CZ)
    , coor_s(CS)
    , py(PY)
    , pz(PZ)
    , weight(W)
    , s_material(std::forward<unique_ptr<Material>>(M)) {}

SectionOS3D::SectionOS3D(const unsigned T, const unsigned MT, const double A, vec&& E)
    : Section(T, SectionType::OS3D, MT, A, std::forward<vec>(E)) {}

/**
 * \brief The deformation is assumed to contain the following.
 *
 *  [0]: u'        \n
 *  [1]: v'        \n
 *  [2]: w'        \n
 *  [3]: v''       \n
 *  [4]: w''       \n
 *  [5]: f         \n
 *  [6]: f'        \n
 *  [7]: f''       \n
 *  [8]: theta_zi  \n
 *  [9]: theta_zj  \n
 *  [10]: theta_yi \n
 *  [11]: theta_yj \n
 *
 */
int SectionOS3D::update_trial_status(const vec& t_deformation) {
    if(const vec incre_deformation = (trial_deformation = t_deformation) - current_deformation; norm(incre_deformation) <= datum::eps) return SUANPAN_SUCCESS;

    const auto& up = trial_deformation(0);
    // const auto& vp = trial_deformation(1);
    // const auto& wp = trial_deformation(2);
    const auto& vpp = trial_deformation(3);
    const auto& wpp = trial_deformation(4);
    const auto& f = trial_deformation(5);
    const auto& fp = trial_deformation(6);
    const auto& fpp = trial_deformation(7);
    const auto theta = trial_deformation.tail(4);

    const rowvec factor = theta.t() * weighing_mat;
    const auto base_strain = dot(factor, theta) + up;

    trial_stiffness.zeros();
    trial_resistance.zeros();
    trial_geometry.zeros();

    for(const auto& I : int_pt) {
        const auto arm_y = I.coor_y - eccentricity(0);
        const auto arm_z = I.coor_z - eccentricity(1);

        if(I.s_material->update_trial_status({base_strain - arm_y * vpp - arm_z * wpp + (arm_z * vpp - arm_y * wpp) * f + .5 * (arm_y * arm_y + arm_z * arm_z) * fp * fp + I.coor_s * fpp, (I.py - arm_z) * fp, (I.pz + arm_y) * fp}) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

        sp_mat de(3, 12);
        de(0, 0) = 1.;
        // de(0, 1) = -arm_y * vp;
        // de(0, 2) = arm_z * wp;
        de(0, 3) = -arm_y + arm_z * f;
        de(0, 4) = -arm_z - arm_y * f;
        de(0, 5) = arm_z * vpp - arm_y * wpp;
        de(0, 6) = (arm_y * arm_y + arm_z * arm_z) * fp;
        de(0, 7) = I.coor_s;
        de(0, 8) = 2. * factor(0);
        de(0, 9) = 2. * factor(1);
        de(0, 10) = 2. * factor(2);
        de(0, 11) = 2. * factor(3);
        de(1, 6) = -arm_z;
        de(2, 6) = arm_y;

        vec shear_scale{1., 1. - I.py / arm_z, 1. + I.pz / arm_y};
        shear_scale(find_nonfinite(shear_scale)).zeros();

        trial_resistance += I.weight * de.t() * I.s_material->get_trial_stress();
        trial_stiffness += I.weight * de.t() * I.s_material->get_trial_stiffness() * diagmat(shear_scale) * de;

        auto axial_force = I.weight * I.s_material->get_trial_stress().at(0);
        const auto major_bending = -arm_y * axial_force, minor_bending = arm_z * axial_force;

        // eq. 7.69 [u',v',w',v'',w'',f,f',f'',theta_zi,theta_zj,theta_yi,theta_yj]
        // trial_geometry(1, 1) += axial_force;
        // trial_geometry(2, 2) += axial_force;
        trial_geometry(3, 5) += minor_bending;
        trial_geometry(5, 3) += minor_bending;
        trial_geometry(4, 5) += major_bending;
        trial_geometry(5, 4) += major_bending;
        trial_geometry(6, 6) += (arm_y * arm_y + arm_z * arm_z) * axial_force;
        axial_force /= 30.;
        trial_geometry(8, 9) -= axial_force;
        trial_geometry(9, 8) -= axial_force;
        trial_geometry(10, 11) -= axial_force;
        trial_geometry(11, 10) -= axial_force;
        axial_force *= 4.;
        trial_geometry(8, 8) += axial_force;
        trial_geometry(9, 9) += axial_force;
        trial_geometry(10, 10) += axial_force;
        trial_geometry(11, 11) += axial_force;
    }

    return SUANPAN_SUCCESS;
}

int SectionOS3D::clear_status() {
    current_deformation = trial_deformation.zeros();
    current_resistance = trial_resistance.zeros();
    current_stiffness = trial_stiffness = initial_stiffness;
    current_geometry = trial_geometry = initial_geometry;
    auto code = 0;
    for(const auto& I : int_pt) code += I.s_material->clear_status();
    return code;
}

int SectionOS3D::commit_status() {
    current_deformation = trial_deformation;
    current_resistance = trial_resistance;
    current_stiffness = trial_stiffness;
    current_geometry = trial_geometry;
    auto code = 0;
    for(const auto& I : int_pt) code += I.s_material->commit_status();
    return code;
}

int SectionOS3D::reset_status() {
    trial_deformation = current_deformation;
    trial_resistance = current_resistance;
    trial_stiffness = current_stiffness;
    trial_geometry = current_geometry;
    auto code = 0;
    for(const auto& I : int_pt) code += I.s_material->reset_status();
    return code;
}

vector<vec> SectionOS3D::record(const OutputType P) {
    if(OutputType::BEAMS == P) {
        vec beam_force(6, fill::zeros);
        for(const auto& I : int_pt) {
            const auto arm_y = I.coor_y - eccentricity(0);
            const auto arm_z = I.coor_z - eccentricity(1);

            const vec force = I.weight * I.s_material->get_current_stress();

            const auto& axial_force = force(0);

            beam_force(0) += axial_force;
            beam_force(1) -= axial_force * arm_y;
            beam_force(2) += axial_force * arm_z;
            beam_force(3) += axial_force * (arm_y * arm_y + arm_z * arm_z);
            beam_force(4) += axial_force * I.coor_s;
            beam_force(5) += arm_y * force(2) - arm_z * force(1);
        }
        return {beam_force};
    }

    return Section::record(P);
}
