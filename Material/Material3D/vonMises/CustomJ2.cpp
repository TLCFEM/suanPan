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

#include "CustomJ2.h"
#include <Domain/DomainBase.h>
#include <Toolbox/utility.h>

double CustomJ2::compute_k(const double p_strain) const { return k_expression->evaluate(p_strain).at(0); }

double CustomJ2::compute_dk(const double p_strain) const { return k_expression->gradient(p_strain).at(0); }

double CustomJ2::compute_h(const double p_strain) const { return h_expression->evaluate(p_strain).at(0); }

double CustomJ2::compute_dh(const double p_strain) const { return h_expression->gradient(p_strain).at(0); }

CustomJ2::CustomJ2(const unsigned T, const double E, const double V, const unsigned K, const unsigned H, const double R)
    : NonlinearJ2(T, E, V, R)
    , k_tag(K)
    , h_tag(H) {}

int CustomJ2::initialize(const shared_ptr<DomainBase>& D) {
    if(!D->find_expression(k_tag)) {
        suanpan_error("Cannot find the assigned expression with tag {}.\n", k_tag);
        return SUANPAN_FAIL;
    }

    k_expression = D->get_expression(k_tag);

    if(k_expression->input_size() != 1 || k_expression->output_size() != 1) {
        suanpan_error("The assigned expression {} is not a univariate function.\n", k_tag);
        return SUANPAN_FAIL;
    }

    if(!D->find_expression(h_tag)) {
        suanpan_error("Cannot find the assigned expression with tag {}.\n", h_tag);
        return SUANPAN_FAIL;
    }

    h_expression = D->get_expression(h_tag);

    if(h_expression->input_size() != 1 || h_expression->output_size() != 1) {
        suanpan_error("The assigned expression {} is not a univariate function.\n", h_tag);
        return SUANPAN_FAIL;
    }

    if(!suanpan::approx_equal(0., h_expression->evaluate(0.).at(0))) {
        suanpan_error("The assigned expression {} does not evaluates to zero.\n", h_tag);
        return SUANPAN_FAIL;
    }

    return NonlinearJ2::initialize(D);
}

unique_ptr<Material> CustomJ2::get_copy() { return make_unique<CustomJ2>(*this); }

void CustomJ2::print() {
    suanpan_info("A 3D J2 hardening model with custom hardening rules.\n");
}
