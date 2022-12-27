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
/**
 * @class B21E
 * @brief The B21E class.
 * @author tlc
 * @date 07/01/2020
 * @version 0.1.0
 * @file B21E.h
 * @addtogroup Beam
 * @ingroup Element
 * @{
 */

#ifndef B21E_H
#define B21E_H

#include "B21.h"

class B21E final : public B21 {
    static const unsigned max_iteration;
    static const double tolerance;

    const uvec a, b;

    vec trial_rotation = zeros(a.n_elem);
    vec current_rotation = zeros(a.n_elem);

public:
    B21E(unsigned,     // tag
         unsigned,     // which
         uvec&&,       // node tags
         unsigned,     // section tag
         unsigned = 6, // integration points
         bool = false  // nonlinear geometry switch
    );

    int update_status() override;

    int commit_status() override;
    int clear_status() override;
    int reset_status() override;
};

#endif

//! @}
