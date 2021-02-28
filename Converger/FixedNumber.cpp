////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2017-2021 Theodore Chang
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
////////////////////////////////////////////////////////////////////////////////

#include "FixedNumber.h"

/**
 * \brief the complete constructor.
 * \param T `unique_tag`
 * \param M `max_itertation`
 * \param P `print_flag`
 */
FixedNumber::FixedNumber(const unsigned T, const unsigned M, const bool P)
	: Converger(T, 1., M, P) {}

/**
 * \brief
 * \return
 */
bool FixedNumber::is_converged() {
	if(is_print()) suanpan_info("iteration counter: %u.\n", ++counter);

	if(get_max_iteration() <= counter) {
		counter = 0;
		set_conv_flag(true);
	}
	else set_conv_flag(false);

	return get_conv_flag();
}
