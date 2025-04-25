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

#include "CustomStrainDegradation.h"

#include <Domain/DomainBase.h>

vec CustomStrainDegradation::compute_positive_degradation(const double t_strain) const { return positive_expression->evaluate(t_strain); }

vec CustomStrainDegradation::compute_negative_degradation(const double t_strain) const {
    if(positive_expression_tag != negative_expression_tag) return negative_expression->evaluate(t_strain);

    auto response = positive_expression->evaluate(fabs(t_strain));
    response(1) = -response(1);
    return response;
}

CustomStrainDegradation::CustomStrainDegradation(const unsigned T, const unsigned MT, const unsigned PET, const unsigned NET)
    : StrainDegradation(T, MT)
    , positive_expression_tag(PET)
    , negative_expression_tag(NET) {}

int CustomStrainDegradation::initialize(const shared_ptr<DomainBase>& D) {
    if(!D->find_expression(positive_expression_tag)) {
        suanpan_error("Cannot find the assigned expression with tag {}.\n", positive_expression_tag);
        return SUANPAN_FAIL;
    }
    if(!D->find_expression(negative_expression_tag)) {
        suanpan_error("Cannot find the assigned expression with tag {}.\n", negative_expression_tag);
        return SUANPAN_FAIL;
    }

    positive_expression = D->get_expression(positive_expression_tag);
    if(positive_expression->input_size() != 1 || positive_expression->output_size() != 2) {
        suanpan_error("The assigned expression {} does not provide two outputs.\n", positive_expression_tag);
        return SUANPAN_FAIL;
    }

    negative_expression = D->get_expression(negative_expression_tag);
    if(negative_expression->input_size() != 1 || negative_expression->output_size() != 2) {
        suanpan_error("The assigned expression {} does not provide two outputs.\n", negative_expression_tag);
        return SUANPAN_FAIL;
    }

    return StrainDegradation::initialize(D);
}

unique_ptr<Material> CustomStrainDegradation::get_copy() { return make_unique<CustomStrainDegradation>(*this); }
