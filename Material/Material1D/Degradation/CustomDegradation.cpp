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

#include "CustomDegradation.h"
#include <Domain/DomainBase.h>

podarray<double> CustomDegradation::compute_degradation(const double t_strain) const { return podarray(expression->evaluate(t_strain).mem, 2); }

CustomDegradation::CustomDegradation(const unsigned T, const unsigned MT, const unsigned ET)
    : Degradation(T, MT)
    , expression_tag(ET) {}

int CustomDegradation::initialize(const shared_ptr<DomainBase>& D) {
    if(!D->find_expression(expression_tag)) {
        suanpan_error("Cannot find the assigned expression with tag {}.\n", expression_tag);
        return SUANPAN_FAIL;
    }

    expression = D->get_expression(expression_tag);

    if(expression->input_size() != 1 || expression->output_size() != 2) {
        suanpan_error("The assigned expression {} does not provide two outputs.\n", expression_tag);
        return SUANPAN_FAIL;
    }

    return Degradation::initialize(D);
}

unique_ptr<Material> CustomDegradation::get_copy() { return make_unique<CustomDegradation>(*this); }
