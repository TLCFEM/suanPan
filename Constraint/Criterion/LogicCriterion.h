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
/**
 * @class LogicCriterion
 * @brief A LogicCriterion class.
 *
 * The LogicCriterion class.
 *
 * @author tlc
 * @date 08/04/2022
 * @version 0.1.0
 * @file LogicCriterion.h
 * @addtogroup Criterion
 * @{
 */

#ifndef LOGICCRITERION_H
#define LOGICCRITERION_H

#include "Criterion.h"

class LogicCriterion : public Criterion {
    const unsigned tag_a, tag_b;

protected:
    shared_ptr<Criterion> criterion_a, criterion_b;

public:
    explicit LogicCriterion(unsigned, unsigned, unsigned, unsigned);

    int initialize(const shared_ptr<DomainBase>&) override;
};

class LogicCriterionAND final : public LogicCriterion {
public:
    using LogicCriterion::LogicCriterion;

    unique_ptr<Criterion> get_copy() override;

    int process(const shared_ptr<DomainBase>&) override;
};

class LogicCriterionOR final : public LogicCriterion {
public:
    using LogicCriterion::LogicCriterion;

    unique_ptr<Criterion> get_copy() override;

    int process(const shared_ptr<DomainBase>&) override;
};

#endif

//! @}
