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

#include "CustomAmplitude.h"
#include <Domain/DomainBase.h>

CustomAmplitude::CustomAmplitude(const unsigned T, const unsigned ET, const unsigned ST)
    : Amplitude(T, ST)
    , e_tag(ET) {}

void CustomAmplitude::initialize(const shared_ptr<DomainBase>& D) {
    if(!D->find<Expression>(e_tag)) {
        suanpan_error("Cannot find expression {}.\n", e_tag);
        D->disable_amplitude(get_tag());
        return;
    }

    expression = D->get<Expression>(e_tag);

    if(expression->input_size() != 1 || expression->output_size() != 1) {
        suanpan_error("The assigned expression {} is not a univariate function.\n", e_tag);
        D->disable_amplitude(get_tag());
    }
}

double CustomAmplitude::get_amplitude(const double T) {
    const auto step_time = T - start_time;

    return step_time <= 0. ? 0. : expression->evaluate(step_time).at(0);
}

void CustomAmplitude::print() {
    suanpan_info("Amplitude using custom function.\n");
}
