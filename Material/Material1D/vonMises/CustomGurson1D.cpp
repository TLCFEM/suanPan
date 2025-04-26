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

#include "CustomGurson1D.h"

#include <Domain/DomainBase.h>

CustomGurson1D::CustomGurson1D(const unsigned T, const unsigned ET, const double E, const double V, const double Q1, const double Q2, const double FN, const double SN, const double EN, const double R)
    : NonlinearGurson1D(T, E, V, Q1, Q2, FN, SN, EN, R)
    , expression_tag(ET) {}

int CustomGurson1D::initialize(const shared_ptr<DomainBase>& D) {
    if(!D->find_expression(expression_tag)) {
        suanpan_error("Cannot find the assigned expression with tag {}.\n", expression_tag);
        return SUANPAN_FAIL;
    }

    expression = D->get_expression(expression_tag);

    if(expression->input_size() != 1 || expression->output_size() != 2) {
        suanpan_error("The assigned expression {} does not generate two outputs.\n", expression_tag);
        return SUANPAN_FAIL;
    }

    return NonlinearGurson1D::initialize(D);
}

vec CustomGurson1D::compute_hardening(const double plastic_strain) const { return expression->evaluate(plastic_strain); }

unique_ptr<Material> CustomGurson1D::get_copy() { return std::make_unique<CustomGurson1D>(*this); }

void CustomGurson1D::print() {
    suanpan_info("A uniaxial Gurson model using custom hardening rule.\n");
    NonlinearGurson1D::print();
}
