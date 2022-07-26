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
/**
 * @class NMB21E
 * @brief The NMB21E class.
 * @author tlc
 * @date 26/07/2022
 * @version 0.1.0
 * @file NMB21E.h
 * @addtogroup Beam
 * @ingroup Element
 * @{
 */

#ifndef NMB21E_H
#define NMB21E_H

#include "NMB21.h"

class NMB21E final : public NMB21 {
    static const unsigned max_iteration;
    static const double tolerance;

    const unsigned which;

    const uvec a{which};
    const uvec b = 1llu == which ? uvec{0, 2} : uvec{0, 1};

    vec trial_rotation = zeros(a.n_elem);
    vec current_rotation = zeros(a.n_elem);

public:
    NMB21E(unsigned,    // tag
           unsigned,    // which
           uvec&&,      // node tags
           unsigned,    // section tag
           bool = false // nonlinear geometry switch
    );

    int update_status() override;

    int commit_status() override;
    int clear_status() override;
    int reset_status() override;
};

#endif

//! @}
