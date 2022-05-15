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
 * @class MaxResistance
 * @brief A MaxResistance class.
 *
 * The MaxResistance class.
 *
 * @author tlc
 * @date 12/02/2020
 * @version 0.1.0
 * @file MaxResistance.h
 * @addtogroup Criterion
 * @{
 */

#ifndef MAXRESISTANCE_H
#define MAXRESISTANCE_H

#include "NodeBasedCriterion.h"

class MaxResistance final : public NodeBasedCriterion {
public:
    using NodeBasedCriterion::NodeBasedCriterion;

    unique_ptr<Criterion> get_copy() override;

    int process(const shared_ptr<DomainBase>&) override;
};

#endif

//! @}
