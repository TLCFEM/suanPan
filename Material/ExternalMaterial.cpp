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

#include "ExternalMaterial.h"

MaterialType ExternalMaterial::get_type(const ExternalMaterialData& D) {
    if(D.size == 6) return MaterialType::D3;
    if(D.size == 3) return MaterialType::D2;
    if(D.size == 1) return MaterialType::D1;
    return MaterialType::D0;
}

ExternalMaterial::ExternalMaterial(const unsigned T, std::vector<double>&& P, void* H)
    : Material(T, MaterialType::D0, 0.)
    , constant(std::move(P))
    // ReSharper disable once CppFunctionalStyleCast
    , cooker(Interface(H)) {
    data.constant = constant.data();
    data.constant_size = static_cast<unsigned>(constant.size());

    int info = ALLOCATE;

    cooker(&data, &info);

    access::rw(density) = data.density;
    access::rw(material_type) = get_type(data);
}

ExternalMaterial::ExternalMaterial(const ExternalMaterial& old_obj)
    : Material(old_obj)
    , constant(old_obj.constant)
    , cooker(old_obj.cooker) {
    data.constant = constant.data();
    data.constant_size = static_cast<unsigned>(constant.size());

    int info = ALLOCATE;

    cooker(&data, &info);

    // need to reinitialize to setup containers
    initialize(nullptr);
}

ExternalMaterial::~ExternalMaterial() {
    int info = DEALLOCATE;

    cooker(&data, &info);
}

bool ExternalMaterial::validate() {
    int info = VALIDATE;

    cooker(&data, &info);

    return 0 == info;
}

int ExternalMaterial::initialize(const shared_ptr<DomainBase>&) {
    PureWrapper(this);

    // ! we need to make all variables wrappers of some external memory

    if(-1 != data.c_strain) current_strain = vec(&data.pool[data.c_strain], data.size, false);
    if(-1 != data.c_strain_rate) current_strain_rate = vec(&data.pool[data.c_strain_rate], data.size, false);
    if(-1 != data.c_stress) current_stress = vec(&data.pool[data.c_stress], data.size, false);

    if(-1 != data.t_strain) trial_strain = vec(&data.pool[data.t_strain], data.size, false);
    if(-1 != data.t_strain_rate) trial_strain_rate = vec(&data.pool[data.t_strain_rate], data.size, false);
    if(-1 != data.t_stress) trial_stress = vec(&data.pool[data.t_stress], data.size, false);

    if(-1 != data.i_stiffness) initial_stiffness = mat(&data.pool[data.i_stiffness], data.size, data.size, false);
    if(-1 != data.c_stiffness) current_stiffness = mat(&data.pool[data.c_stiffness], data.size, data.size, false);
    if(-1 != data.t_stiffness) trial_stiffness = mat(&data.pool[data.t_stiffness], data.size, data.size, false);

    if(-1 != data.i_damping) initial_damping = mat(&data.pool[data.i_damping], data.size, data.size, false);
    if(-1 != data.c_damping) current_damping = mat(&data.pool[data.c_damping], data.size, data.size, false);
    if(-1 != data.t_damping) trial_damping = mat(&data.pool[data.t_damping], data.size, data.size, false);

    return SUANPAN_SUCCESS;
}

void ExternalMaterial::initialize_history(unsigned) {}

void ExternalMaterial::set_initial_history(const vec&) {}

unique_ptr<Material> ExternalMaterial::get_copy() { return make_unique<ExternalMaterial>(*this); }

int ExternalMaterial::update_trial_status(const vec& t_strain) {
    if(-1 == data.t_strain) return SUANPAN_FAIL;

    trial_strain = t_strain;

    int info = STATIC_UPDATE;

    cooker(&data, &info);

    return info;
}

int ExternalMaterial::update_trial_status(const vec& t_strain, const vec& t_strain_rate) {
    if(-1 == data.t_strain || -1 == data.t_strain_rate) return SUANPAN_FAIL;

    trial_strain = t_strain;
    trial_strain_rate = t_strain_rate;

    int info = DYNAMIC_UPDATE;

    cooker(&data, &info);

    return info;
}

int ExternalMaterial::commit_status() {
    int info = COMMIT;

    cooker(&data, &info);

    return info;
}

int ExternalMaterial::reset_status() {
    int info = RESET;

    cooker(&data, &info);

    return info;
}

int ExternalMaterial::clear_status() {
    int info = CLEAR;

    cooker(&data, &info);

    return info;
}

std::vector<vec> ExternalMaterial::record(OutputType) { return {}; }
