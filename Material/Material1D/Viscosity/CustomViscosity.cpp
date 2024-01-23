/*******************************************************************************
 * Copyright (C) 2017-2024 Theodore Chang
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

#include "CustomViscosity.h"
#include <Domain/DomainBase.h>

double CustomViscosity::compute_du(const double strain, const double strain_rate) const { return expression->evaluate(vec{strain, strain_rate}).at(1); }

double CustomViscosity::compute_dv(const double strain, const double strain_rate) const { return expression->evaluate(vec{strain, strain_rate}).at(2); }

double CustomViscosity::compute_damping_coefficient(const double strain, const double strain_rate) const { return expression->evaluate(vec{strain, strain_rate}).at(0); }

CustomViscosity::CustomViscosity(const unsigned T, const unsigned ET)
    : NonlinearViscosity(T, 0., 1.)
    , expression_tag(ET) {}

int CustomViscosity::initialize(const shared_ptr<DomainBase>& D) {
    if(!D->find_expression(expression_tag)) {
        suanpan_error("Cannot find the assigned expression with tag {}.\n", expression_tag);
        return SUANPAN_FAIL;
    }

    expression = D->get_expression(expression_tag);

    if(expression->input_size() != 2 || expression->output_size() != 3) {
        suanpan_error("An expression with two inputs and three outputs is required.\n");
        return SUANPAN_FAIL;
    }

    return NonlinearViscosity::initialize(D);
}

unique_ptr<Material> CustomViscosity::get_copy() { return make_unique<CustomViscosity>(*this); }
