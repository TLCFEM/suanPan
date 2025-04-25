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

#include "CustomMises1D.h"

#include <Domain/DomainBase.h>
#include <Toolbox/utility.h>

double CustomMises1D::compute_k(const double p_strain) const { return k_expression->evaluate(p_strain).at(0); }

double CustomMises1D::compute_dk(const double p_strain) const { return k_expression->gradient(p_strain).at(0); }

double CustomMises1D::compute_h(const double p_strain) const { return h_expression->evaluate(p_strain).at(0); }

double CustomMises1D::compute_dh(const double p_strain) const { return h_expression->gradient(p_strain).at(0); }

CustomMises1D::CustomMises1D(const unsigned T, const double E, const unsigned K, const unsigned H, const double R)
    : NonlinearMises1D(T, E, R)
    , k_tag(K)
    , h_tag(H) {}

int CustomMises1D::initialize(const shared_ptr<DomainBase>& D) {
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

    return NonlinearMises1D::initialize(D);
}

unique_ptr<Material> CustomMises1D::get_copy() { return make_unique<CustomMises1D>(*this); }
