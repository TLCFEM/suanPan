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

#include "CustomHoffman.h"

#include <Domain/DomainBase.h>
#include <Toolbox/utility.h>

double CustomHoffman::compute_k(const double p_strain) const { return k_expression->evaluate(p_strain).at(0); }

double CustomHoffman::compute_dk(const double p_strain) const { return k_expression->gradient(p_strain).at(0); }

CustomHoffman::CustomHoffman(const unsigned T, vec&& E, vec&& V, vec&& S, const unsigned KT, const double R)
    : NonlinearOrthotropic(T, OrthotropicType::Hoffman, std::move(E), std::move(V), std::move(S), R)
    , k_tag(KT) {}

int CustomHoffman::initialize(const shared_ptr<DomainBase>& D) {
    if(!D->find_expression(k_tag)) {
        suanpan_error("Cannot find the assigned expression with tag {}.\n", k_tag);
        return SUANPAN_FAIL;
    }

    k_expression = D->get_expression(k_tag);

    if(k_expression->input_size() != 1 || k_expression->output_size() != 1) {
        suanpan_error("The assigned expression {} is not a univariate function.\n", k_tag);
        return SUANPAN_FAIL;
    }

    if(!suanpan::approx_equal(1., k_expression->evaluate(0.).at(0))) {
        suanpan_error("The assigned expression {} does not evaluates to unity for trial plastic strain.\n", k_tag);
        return SUANPAN_FAIL;
    }

    return NonlinearOrthotropic::initialize(D);
}

unique_ptr<Material> CustomHoffman::get_copy() { return std::make_unique<CustomHoffman>(*this); }

void CustomHoffman::print() { suanpan_info("A 3D nonlinear Hoffman model using custom hardening function.\n"); }
