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

#include "Homogeneous.h"

#include <Domain/DomainBase.h>
#include <Material/Material2D/Material2D.h>
#include <Toolbox/IntegrationPlan.h>

Homogeneous::IntegrationPoint::IntegrationPoint(const double E, const double F, unique_ptr<Material>&& M)
    : eccentricity(E)
    , factor(F)
    , s_material(std::move(M)) {}

Homogeneous::IntegrationPoint::IntegrationPoint(const IntegrationPoint& old_obj)
    : eccentricity(old_obj.eccentricity)
    , factor(old_obj.factor)
    , s_material(suanpan::make_copy(old_obj.s_material)) {}

Homogeneous::Homogeneous(const unsigned T, const unsigned MT, const double TH, const unsigned IP)
    : SectionShell(T, MT)
    , num_ip(IP > 20 ? 20 : IP)
    , thickness(TH) {}

int Homogeneous::initialize(const shared_ptr<DomainBase>& D) {
    const auto mat_proto = std::dynamic_pointer_cast<Material2D>(D->get_material(material_tag));

    if(nullptr == mat_proto) return SUANPAN_FAIL;

    const IntegrationPlan plan(1, num_ip, IntegrationType::GAUSS);

    initial_membrane_stiffness.zeros(3, 3);
    initial_plate_stiffness.zeros(3, 3);

    for(unsigned I = 0; I < plan.n_rows; ++I) {
        int_pt.emplace_back(.5 * thickness * plan(I, 0), plan(I, 1), mat_proto->get_copy());

        const auto& c_pt = int_pt.back();

        initial_membrane_stiffness += c_pt.factor * c_pt.s_material->get_trial_stiffness();
        initial_plate_stiffness += c_pt.factor * c_pt.eccentricity * c_pt.eccentricity * c_pt.s_material->get_trial_stiffness();
    }

    trial_membrane_stiffness = current_membrane_stiffness = initial_membrane_stiffness;
    trial_plate_stiffness = current_membrane_stiffness = initial_membrane_stiffness;

    return SUANPAN_SUCCESS;
}

unique_ptr<SectionShell> Homogeneous::get_copy() { return make_unique<Homogeneous>(*this); }

int Homogeneous::update_trial_status(const vec& m_strain, const vec& p_strain) {
    trial_membrane_strain = m_strain;
    trial_plate_strain = p_strain;

    trial_membrane_stress.zeros(3);
    trial_plate_stress.zeros(3);
    trial_membrane_stiffness.zeros(3, 3);
    trial_plate_stiffness.zeros(3, 3);

    for(const auto& I : int_pt) {
        if(I.s_material->update_trial_status(m_strain + I.eccentricity * p_strain) != SUANPAN_SUCCESS) return SUANPAN_FAIL;
        trial_membrane_stress += I.factor * I.s_material->get_trial_stress();
        trial_membrane_stiffness += I.factor * I.s_material->get_trial_stiffness();
        trial_plate_stress += I.factor * I.eccentricity * I.s_material->get_trial_stress();
        trial_plate_stiffness += I.factor * I.eccentricity * I.eccentricity * I.s_material->get_trial_stiffness();
    }

    return SUANPAN_SUCCESS;
}

int Homogeneous::clear_status() {
    current_membrane_strain.zeros();
    current_membrane_strain_rate.zeros();
    current_membrane_stress.zeros();
    current_membrane_stiffness = initial_membrane_stiffness;
    current_plate_strain.zeros();
    current_plate_strain_rate.zeros();
    current_plate_stress.zeros();
    current_plate_stiffness = initial_plate_stiffness;
    auto code = reset_status();
    for(const auto& I : int_pt) code += I.s_material->clear_status();
    return code;
}

int Homogeneous::commit_status() {
    current_membrane_strain = trial_membrane_strain;
    current_membrane_strain_rate = trial_membrane_strain_rate;
    current_membrane_stress = trial_membrane_stress;
    current_membrane_stiffness = trial_membrane_stiffness;
    current_plate_strain = trial_plate_strain;
    current_plate_strain_rate = trial_plate_strain_rate;
    current_plate_stress = trial_plate_stress;
    current_plate_stiffness = trial_plate_stiffness;
    auto code = 0;
    for(const auto& I : int_pt) code += I.s_material->commit_status();
    return code;
}

int Homogeneous::reset_status() {
    trial_membrane_strain = current_membrane_strain;
    trial_membrane_strain_rate = current_membrane_strain_rate;
    trial_membrane_stress = current_membrane_stress;
    trial_membrane_stiffness = current_membrane_stiffness;
    trial_plate_strain = current_plate_strain;
    trial_plate_strain_rate = current_plate_strain_rate;
    trial_plate_stress = current_plate_stress;
    trial_plate_stiffness = current_plate_stiffness;
    auto code = 0;
    for(const auto& I : int_pt) code += I.s_material->reset_status();
    return code;
}
