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

#include "AbsResidual.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>

AbsResidual::AbsResidual(const unsigned T, const double E, const unsigned M, const bool P)
    : Converger(T, E, M, P) {}

unique_ptr<Converger> AbsResidual::get_copy() { return make_unique<AbsResidual>(*this); }

bool AbsResidual::is_converged(unsigned) {
    auto& W = get_domain().lock()->get_factory();

    set_error(norm(get_residual()) / static_cast<double>(W->get_size()));
    set_conv_flag(get_tolerance() > get_error());

    if(is_print()) sp_info("--> Absolute Residual: {:.5E}.\n", get_error());

    return get_conv_flag();
}
