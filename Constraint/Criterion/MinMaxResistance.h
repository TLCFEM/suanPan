/*******************************************************************************
 * Copyright (C) 2017-2026 Theodore Chang
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
 * @class ResistanceCriterion
 * @brief A ResistanceCriterion class.
 *
 * The ResistanceCriterion class.
 *
 * @author tlc
 * @date 29/07/2026
 * @version 1.0.0
 * @file MinMaxResistance.h
 * @addtogroup Criterion
 * @{
 */

#ifndef MINMAXRESISTANCE_H
#define MINMAXRESISTANCE_H

#include "NodeBasedCriterion.h"

class ResistanceCriterion : public NodeBasedCriterion {
    [[nodiscard]] virtual bool check(double) const = 0;

public:
    using NodeBasedCriterion::NodeBasedCriterion;

    int process(const shared_ptr<DomainBase>&) final;
};

class MaxResistance final : public ResistanceCriterion {
    [[nodiscard]] bool check(double) const override;

public:
    using ResistanceCriterion::ResistanceCriterion;

    unique_ptr<Criterion> unique_copy() override;
};

class MinResistance final : public ResistanceCriterion {
    [[nodiscard]] bool check(double) const override;

public:
    using ResistanceCriterion::ResistanceCriterion;

    unique_ptr<Criterion> unique_copy() override;
};

#endif

//! @}
