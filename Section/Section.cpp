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

#include "Section.h"
#include <Domain/DomainBase.h>
#include <Material/Material.h>
#include <Recorder/OutputType.h>

Section::Section(const unsigned T, const SectionType ST, const unsigned MT, const double A, vec&& EC)
    : DataSection{MT, ST, EC.head(2), A}
    , CopiableTag(T) {}

SectionType Section::get_section_type() const { return section_type; }

double Section::get_area() const { return area; }

double Section::get_linear_density() const { return linear_density; }

int Section::initialize_base(const shared_ptr<DomainBase>& D) {
    if(initialized) return SUANPAN_SUCCESS;

    if(0u != material_tag) {
        if(!D->find<Material>(material_tag)) {
            suanpan_warning("Section {} disabled as material {} cannot be found.\n", get_tag(), material_tag);
            return SUANPAN_FAIL;
        }
        if(SectionType::OS3D == section_type) {
            if(MaterialType::OS != D->get<Material>(material_tag)->get_material_type()) {
                suanpan_warning("Section {} disabled as material {} has a wrong type, use OS type material only.\n", get_tag(), material_tag);
                return SUANPAN_FAIL;
            }
        }
        else if(MaterialType::D1 != D->get<Material>(material_tag)->get_material_type()) {
            suanpan_warning("Section {} disabled as material {} has a wrong type, use 1D material only.\n", get_tag(), material_tag);
            return SUANPAN_FAIL;
        }
    }

    const auto size = static_cast<unsigned>(section_type);

    if(current_deformation.is_empty()) current_deformation.zeros(size);
    if(trial_deformation.is_empty()) trial_deformation.zeros(size);

    // current_deformation_rate.zeros(size);
    // trial_deformation_rate.zeros(size);

    if(current_resistance.is_empty()) current_resistance.zeros(size);
    if(trial_resistance.is_empty()) trial_resistance.zeros(size);

    if(initial_stiffness.is_empty()) initial_stiffness.zeros(size, size);
    if(trial_stiffness.is_empty()) trial_stiffness.zeros(size, size);
    if(current_stiffness.is_empty()) current_stiffness.zeros(size, size);

    return SUANPAN_SUCCESS;
}

void Section::set_initialized(const bool F) const { access::rw(initialized) = F; }

void Section::set_symmetric(const bool F) const { access::rw(symmetric) = F; }

bool Section::is_initialized() const { return initialized; }

bool Section::is_symmetric() const { return symmetric; }

void Section::set_eccentricity(const vec& E) const { access::rw(eccentricity) = E; }

const vec& Section::get_eccentricity() const { return eccentricity; }

void Section::set_characteristic_length(const double L) const { access::rw(characteristic_length) = std::max(datum::eps, L); }

double Section::get_characteristic_length() const { return characteristic_length; }

const vec& Section::get_trial_deformation() const { return trial_deformation; }

const vec& Section::get_trial_deformation_rate() const { return trial_deformation_rate; }

const vec& Section::get_trial_resistance() const { return trial_resistance; }

const mat& Section::get_trial_stiffness() const { return trial_stiffness; }

const mat& Section::get_trial_geometry() const { return trial_geometry; }

const vec& Section::get_current_deformation() const { return current_deformation; }

const vec& Section::get_current_deformation_rate() const { return current_deformation_rate; }

const vec& Section::get_current_resistance() const { return current_resistance; }

const mat& Section::get_current_stiffness() const { return current_stiffness; }

const mat& Section::get_current_geometry() const { return current_geometry; }

const mat& Section::get_initial_stiffness() const { return initial_stiffness; }

const mat& Section::get_initial_geometry() const { return initial_geometry; }

int Section::update_incre_status(const double i_strain) {
    const vec i_vec_strain{i_strain};
    return update_incre_status(i_vec_strain);
}

int Section::update_incre_status(const double i_strain, const double i_strain_rate) {
    const vec i_vec_strain{i_strain};
    const vec i_vec_strain_rate{i_strain_rate};
    return update_incre_status(i_vec_strain, i_vec_strain_rate);
}

int Section::update_trial_status(const double t_strain) {
    const vec t_vec_strain{t_strain};
    return update_trial_status(t_vec_strain);
}

int Section::update_trial_status(const double t_strain, const double t_strain_rate) {
    const vec t_vec_strain{t_strain};
    const vec t_vec_strain_rate{t_strain_rate};
    return update_trial_status(t_vec_strain, t_vec_strain_rate);
}

int Section::update_incre_status(const vec& i_deformation) { return update_trial_status(current_deformation + i_deformation); }

int Section::update_incre_status(const vec& i_deformation, const vec& i_deformation_rate) { return update_trial_status(current_deformation + i_deformation, current_deformation_rate + i_deformation_rate); }

int Section::update_trial_status(const vec&) { throw invalid_argument("hidden method called"); }

int Section::update_trial_status(const vec& t_deformation, const vec&) { return update_trial_status(t_deformation); }

std::vector<vec> Section::record(const OutputType P) {
    if(P == OutputType::E) return {current_deformation};
    if(P == OutputType::S) return {current_resistance};
    if(P == OutputType::PE) return {current_deformation - solve(initial_stiffness, current_resistance)};
    if(P == OutputType::EE) return {solve(initial_stiffness, current_resistance)};

    return {};
}

unique_ptr<Section> suanpan::make_copy(const shared_ptr<Section>& S) { return S->get_copy(); }

unique_ptr<Section> suanpan::make_copy(const unique_ptr<Section>& S) { return S->get_copy(); }

unique_ptr<Section> suanpan::initialized_section_copy(const shared_ptr<DomainBase>& D, const uword T) {
    if(!D->find<Section>(T)) return nullptr;

    auto copy = D->get<Section>(T)->get_copy();

    if(copy->is_initialized()) return copy;

    if(SUANPAN_SUCCESS != copy->initialize_base(D) || SUANPAN_SUCCESS != copy->initialize(D)) return nullptr;

    copy->set_initialized(true);

    return copy;
}
