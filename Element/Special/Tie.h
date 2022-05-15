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
 * @class Tie
 * @brief A Tie class.
 *
 * The Tie class.
 *
 * @author tlc
 * @date 29/07/2017
 * @version 0.1.0
 * @file Tie.h
 * @addtogroup Constraint
 * @{
 */

#ifndef TIE_H
#define TIE_H

#include <Element/Element.h>

class Tie final : public Element {
    const uvec dof_pool;
    const vec weight_pool;
    const double pseudo_load;
    const double penalty;
    uvec t_span;

public:
    Tie(unsigned, uvec&&, uvec&&, vec&&, double, double = 1E8);

    int initialize(const shared_ptr<DomainBase>&) override;

    int update_status() override;

    int clear_status() override;
    int commit_status() override;
    int reset_status() override;
};

#endif

//! @}
