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

#include "CustomElastic1D.h"
#include <Domain/DomainBase.h>

CustomElastic1D::CustomElastic1D(const unsigned T, const unsigned ET, const double R)
    : Material1D(T, R)
    , expression_tag(ET) {}

int CustomElastic1D::initialize(const shared_ptr<DomainBase>& D) {
    if(!D->find_expression(expression_tag)) {
        suanpan_error("Cannot find the assigned expression with tag {}.\n", expression_tag);
        return SUANPAN_FAIL;
    }

    expression = D->get_expression(expression_tag);

    if(expression->size() != 1) {
        suanpan_error("The assigned expression {} is not a univariate function.\n", expression_tag);
        return SUANPAN_FAIL;
    }

    trial_stiffness = current_stiffness = initial_stiffness = expression->gradient(vec{0.}).at(0);

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> CustomElastic1D::get_copy() { return make_unique<CustomElastic1D>(*this); }

int CustomElastic1D::update_trial_status(const vec& t_strain) {
    trial_stress = expression->evaluate(trial_strain = t_strain);
    trial_stiffness = expression->gradient(trial_strain).at(0);
    return SUANPAN_SUCCESS;
}

int CustomElastic1D::clear_status() {
    current_strain = trial_strain.zeros();
    current_stress = trial_stress.zeros();
    current_stiffness = trial_stiffness = initial_stiffness;
    return 0;
}

int CustomElastic1D::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_stiffness = trial_stiffness;
    return 0;
}

int CustomElastic1D::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_stiffness = current_stiffness;
    return 0;
}

void CustomElastic1D::print() {
    suanpan_info("A uniaxial elastic model using custom constitutive equation.\n");
    if(expression) expression->print();
}
