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

#include "PlaneStrain.h"
#include <Domain/DomainBase.h>

const uvec PlaneStrain::F{0, 1, 3};

PlaneStrain::PlaneStrain(const unsigned T, const unsigned BT, const unsigned ST)
    : Material2D(T, PlaneType::E, 0.)
    , FA(0 == ST ? std::initializer_list<uword>{} : 1 == ST ? std::initializer_list<uword>{0, 2} : std::initializer_list<uword>{1, 2})
    , FB(0 == ST ? std::initializer_list<uword>{} : 1 == ST ? std::initializer_list<uword>{2, 5} : std::initializer_list<uword>{2, 4})
    , base_tag(BT) {}

int PlaneStrain::initialize(const shared_ptr<DomainBase>& D) {
    base = suanpan::initialized_material_copy(D, base_tag);

    if(nullptr == base || base->get_material_type() != MaterialType::D3) {
        suanpan_error("A valid 3D host material is required.\n");
        return SUANPAN_FAIL;
    }

    access::rw(density) = base->get_density();

    current_stiffness = trial_stiffness = initial_stiffness = base->get_initial_stiffness()(F, F);

    return SUANPAN_SUCCESS;
}

double PlaneStrain::get_parameter(const ParameterType P) const { return base->get_parameter(P); }

unique_ptr<Material> PlaneStrain::get_copy() { return make_unique<PlaneStrain>(*this); }

int PlaneStrain::update_trial_status(const vec& t_strain) {
    vec full_strain(6, fill::zeros);

    full_strain(F) = trial_strain = t_strain;
    if(!FA.empty()) full_strain(FB) = trial_strain(FA);

    if(SUANPAN_SUCCESS != base->update_trial_status(full_strain)) return SUANPAN_FAIL;

    trial_stress = base->get_trial_stress()(F);

    trial_stiffness = base->get_trial_stiffness()(F, F);

    return SUANPAN_SUCCESS;
}

int PlaneStrain::clear_status() {
    current_strain = trial_strain.zeros();
    current_stress = trial_stress.zeros();
    current_stiffness = trial_stiffness = initial_stiffness;
    return base->clear_status();
}

int PlaneStrain::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_stiffness = trial_stiffness;
    return base->commit_status();
}

int PlaneStrain::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_stiffness = current_stiffness;
    return base->reset_status();
}

vector<vec> PlaneStrain::record(const OutputType P) { return base->record(P); }

void PlaneStrain::print() {
    suanpan_info("A plane strain wrapper.\n");
    if(base) base->print();
}
