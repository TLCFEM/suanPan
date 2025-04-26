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

#include "CustomDP.h"

#include <Domain/DomainBase.h>

double CustomDP::compute_c(const double p_strain) const { return c_expression->evaluate(p_strain).at(0); }

double CustomDP::compute_dc(const double p_strain) const { return c_expression->gradient(p_strain).at(0); }

CustomDP::CustomDP(const unsigned T, const double E, const double V, const double ETAY, const double ETAF, const double XI, const unsigned CT, const double R)
    : NonlinearDruckerPrager(T, E, V, ETAY, ETAF, XI, R)
    , c_tag(CT) {}

int CustomDP::initialize(const shared_ptr<DomainBase>& D) {
    if(!D->find_expression(c_tag)) {
        suanpan_error("Cannot find the assigned expression with tag {}.\n", c_tag);
        return SUANPAN_FAIL;
    }

    c_expression = D->get_expression(c_tag);

    if(c_expression->input_size() != 1 || c_expression->output_size() != 1) {
        suanpan_error("The assigned expression {} is not a univariate function.\n", c_tag);
        return SUANPAN_FAIL;
    }

    return NonlinearDruckerPrager::initialize(D);
}

unique_ptr<Material> CustomDP::get_copy() { return std::make_unique<CustomDP>(*this); }

void CustomDP::print() {
    suanpan_info("A 3D nonlinear model using Drucker-Prager yielding criterion with custom cohesion function.\n");
}
