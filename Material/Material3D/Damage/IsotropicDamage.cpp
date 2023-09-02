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

#include "IsotropicDamage.h"
#include <Domain/DomainBase.h>

IsotropicDamage::IsotropicDamage(const unsigned T, const unsigned MT)
	: Material3D(T, 0.)
	, mat_tag(MT) {}

IsotropicDamage::IsotropicDamage(const IsotropicDamage& old_obj)
	: Material3D(old_obj)
	, mat_tag(old_obj.mat_tag)
	, mat_ptr(suanpan::make_copy(old_obj.mat_ptr)) {}

int IsotropicDamage::initialize(const shared_ptr<DomainBase>& D) {
	mat_ptr = suanpan::initialized_material_copy(D, mat_tag);

	if(nullptr == mat_ptr || mat_ptr->get_material_type() != MaterialType::D3) return SUANPAN_FAIL;

	access::rw(density) = mat_ptr->get_parameter(ParameterType::DENSITY);

	trial_stiffness = current_stiffness = initial_stiffness = mat_ptr->get_initial_stiffness();

	return SUANPAN_SUCCESS;
}

double IsotropicDamage::get_parameter(const ParameterType P) const { return mat_ptr->get_parameter(P); }

int IsotropicDamage::update_trial_status(const vec& t_strain) {
	if(mat_ptr->update_trial_status(trial_strain = t_strain) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

	trial_stress = mat_ptr->get_trial_stress();
	trial_stiffness = mat_ptr->get_trial_stiffness();

	compute_damage();

	return SUANPAN_SUCCESS;
}

int IsotropicDamage::clear_status() {
	trial_strain = current_strain.zeros();
	trial_stress = current_stress.zeros();
	trial_history = current_history = initial_history;
	trial_stiffness = current_stiffness = initial_stiffness;

	return mat_ptr->clear_status();
}

int IsotropicDamage::commit_status() {
	current_strain = trial_strain;
	current_stress = trial_stress;
	current_history = trial_history;
	current_stiffness = trial_stiffness;

	return mat_ptr->commit_status();
}

int IsotropicDamage::reset_status() {
	trial_strain = current_strain;
	trial_stress = current_stress;
	trial_history = current_history;
	trial_stiffness = current_stiffness;

	return mat_ptr->reset_status();
}
