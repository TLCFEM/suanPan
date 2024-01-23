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

#include "RelIncreEnergy.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>

/**
 * \brief The complete constructor.
 * \param T `unique_tag`
 * \param E `tolerance`
 * \param M `max_iteration`
 * \param P `print_flag`
 */
RelIncreEnergy::RelIncreEnergy(const unsigned T, const double E, const unsigned M, const bool P)
    : Converger(T, E, M, P) {}

unique_ptr<Converger> RelIncreEnergy::get_copy() { return make_unique<RelIncreEnergy>(*this); }

bool RelIncreEnergy::is_converged(const unsigned counter) {
    const auto D = get_domain().lock();
    auto& W = D->get_factory();

    const auto energy = fabs(dot(W->get_ninja(), get_residual()));
    if(0u == counter) ref_energy = energy;
    set_error(energy / ref_energy);
    set_conv_flag(get_tolerance() > get_error());

    if(is_print())
        suanpan_info("-> Relative Incremental Energy: {:.5E}.\n", get_error());

    return get_conv_flag();
}
