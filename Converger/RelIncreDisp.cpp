/*******************************************************************************
 * Copyright (C) 2017-2022 Theodore Chang
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

#include "RelIncreDisp.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>

RelIncreDisp::RelIncreDisp(const unsigned T, const double E, const unsigned M, const bool P)
    : Converger(T, E, M, P) {}

unique_ptr<Converger> RelIncreDisp::get_copy() { return make_unique<RelIncreDisp>(*this); }

bool RelIncreDisp::is_converged(unsigned) {
    auto& W = get_domain().lock()->get_factory();

    const auto rel_incre_disp = norm(W->get_ninja()) / norm(W->get_incre_displacement() + W->get_ninja());
    set_error(std::isfinite(rel_incre_disp) ? rel_incre_disp : 1.);
    set_conv_flag(get_tolerance() > get_error());

    if(is_print()) suanpan_info("relative incremental displacement error: %.5E.\n", get_error());

    return get_conv_flag();
}
