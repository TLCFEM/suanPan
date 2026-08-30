/*******************************************************************************
 * Copyright (C) 2017-2026 Theodore Chang
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

    if(expression->input_size() != 1 || expression->output_size() != 1) {
        suanpan_error("The assigned expression {} is not a univariate function.\n", expression_tag);
        return SUANPAN_FAIL;
    }

    trial_stiffness = current_stiffness = initial_stiffness = expression->gradient(vec{0.});

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> CustomElastic1D::unique_copy() { return std::make_unique<CustomElastic1D>(*this); }

int CustomElastic1D::update_trial_status(const vec& t_strain) {
    trial_stress = expression->evaluate(trial_strain = t_strain);
    trial_stiffness = expression->gradient(trial_strain);
    return SUANPAN_SUCCESS;
}

void CustomElastic1D::print() {
    suanpan_info("A uniaxial elastic model using custom constitutive equation.\n");
    if(expression) expression->print();
}
