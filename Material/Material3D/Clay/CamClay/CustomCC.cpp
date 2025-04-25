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

#include "CustomCC.h"

#include <Domain/DomainBase.h>

double CustomCC::compute_a(const double hardening) const { return a_expression->evaluate(hardening).at(0); }

double CustomCC::compute_da(const double hardening) const { return a_expression->gradient(hardening).at(0); }

CustomCC::CustomCC(const unsigned T, const double E, const double V, const double B, const double M, const double P, const unsigned AT, const double R)
    : NonlinearCamClay(T, E, V, B, M, P, R)
    , a_tag(AT) {}

int CustomCC::initialize(const shared_ptr<DomainBase>& D) {
    if(!D->find_expression(a_tag)) {
        suanpan_error("Cannot find the assigned expression with tag {}.\n", a_tag);
        return SUANPAN_FAIL;
    }

    a_expression = D->get_expression(a_tag);

    if(a_expression->input_size() != 1 || a_expression->output_size() != 1) {
        suanpan_error("The assigned expression {} is not a univariate function.\n", a_tag);
        return SUANPAN_FAIL;
    }

    return NonlinearCamClay::initialize(D);
}

unique_ptr<Material> CustomCC::get_copy() { return make_unique<CustomCC>(*this); }

void CustomCC::print() {
    suanpan_info("A 3D Cam-Clay model using custom hardening.\n");
}
