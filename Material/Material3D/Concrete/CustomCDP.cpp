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

#include "CustomCDP.h"
#include <Domain/DomainBase.h>
#include <Toolbox/utility.h>

pod6 CustomCDP::compute_tension_backbone(const double kappa) const {
    pod6 response;
    const auto t_response = t_expression->evaluate(kappa);
    for(auto I = 0llu; I < t_response.n_elem; ++I) response[I] = t_response(I);

    if(response[1] < 0.) {
        response[1] = -response[1];
        response[4] = -response[4];
    }
    if(response[2] < 0.) {
        response[2] = -response[2];
        response[5] = -response[5];
    }
    return response;
}

pod6 CustomCDP::compute_compression_backbone(const double kappa) const {
    pod6 response;
    const auto c_response = c_expression->evaluate(kappa);
    for(auto I = 0llu; I < c_response.n_elem; ++I) response[I] = c_response(I);

    if(response[1] > 0.) {
        response[1] = -response[1];
        response[4] = -response[4];
    }
    if(response[2] > 0.) {
        response[2] = -response[2];
        response[5] = -response[5];
    }
    return response;
}

CustomCDP::CustomCDP(const unsigned T, const unsigned TT, const unsigned CT, const double E, const double V, const double GT, const double GC, const double AP, const double BC, const double S, const double R)
    : NonlinearCDP(T, E, V, GT, GC, AP, BC, S, R)
    , t_tag(TT)
    , c_tag(CT) {}

int CustomCDP::initialize(const shared_ptr<DomainBase>& D) {
    if(!D->find_expression(t_tag)) {
        suanpan_error("Cannot find the assigned expression with tag {}.\n", t_tag);
        return SUANPAN_FAIL;
    }

    t_expression = D->get_expression(t_tag);
    if(t_expression->input_size() != 1 || t_expression->output_size() != 6) {
        suanpan_error("The assigned expression {} does not meet requirements: 1 input, 6 outputs.\n", t_tag);
        return SUANPAN_FAIL;
    }

    if(!D->find_expression(c_tag)) {
        suanpan_error("Cannot find the assigned expression with tag {}.\n", c_tag);
        return SUANPAN_FAIL;
    }

    c_expression = D->get_expression(c_tag);
    if(c_expression->input_size() != 1 || c_expression->output_size() != 6) {
        suanpan_error("The assigned expression {} does not meet requirements: 1 input, 6 outputs.\n", c_tag);
        return SUANPAN_FAIL;
    }

    return NonlinearCDP::initialize(D);
}

unique_ptr<Material> CustomCDP::get_copy() { return make_unique<CustomCDP>(*this); }
