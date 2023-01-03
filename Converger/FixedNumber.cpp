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

#include "FixedNumber.h"

/**
 * \brief the complete constructor.
 * \param T `unique_tag`
 * \param M `max_iteration`
 * \param P `print_flag`
 */
FixedNumber::FixedNumber(const unsigned T, const unsigned M, const bool P)
    : Converger(T, 1., M, P) {}

unique_ptr<Converger> FixedNumber::get_copy() { return make_unique<FixedNumber>(*this); }

bool FixedNumber::is_converged(const unsigned counter) {
    if(is_print())
        suanpan_info("--> Iteration Counter: {}.\n", counter);

    set_conv_flag(get_max_iteration() <= counter);

    return get_conv_flag();
}
