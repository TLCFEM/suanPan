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

#include "AbsIncreEnergy.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>

/**
 * \brief The complete constructor.
 * \param T `unique_tag`
 * \param E `tolerance`
 * \param M `max_iteration`
 * \param P `print_flag`
 */
AbsIncreEnergy::AbsIncreEnergy(const unsigned T, const double E, const unsigned M, const bool P)
    : Converger(T, E, M, P) {}

unique_ptr<Converger> AbsIncreEnergy::get_copy() { return make_unique<AbsIncreEnergy>(*this); }

bool AbsIncreEnergy::is_converged() {
    const auto& D = get_domain().lock();

    if(auto& W = D->get_factory(); W->get_reference_load().is_empty() || W->get_trial_load_factor().is_empty()) set_error(fabs(dot(W->get_ninja(), W->get_trial_load() - W->get_sushi())) / static_cast<double>(W->get_ninja().n_elem));
    else set_error(fabs(dot(W->get_ninja(), W->get_reference_load() * W->get_trial_load_factor() + W->get_trial_load() - W->get_sushi())) / static_cast<double>(W->get_ninja().n_elem));

    set_conv_flag(get_tolerance() > get_error());

    if(is_print()) suanpan_info("absolute energy increment error: %.5E.\n", get_error());

    return get_conv_flag();
}
