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

    [[nodiscard]] virtual int check(int, int) const = 0;

protected:
    shared_ptr<Criterion> criterion_a, criterion_b;

public:
    LogicCriterion(unsigned, unsigned, unsigned);

    int initialize(const shared_ptr<DomainBase>&) override;

    int process(const shared_ptr<DomainBase>&) final;
};

class LogicCriterionAND final : public LogicCriterion {
    [[nodiscard]] int check(int, int) const override;

public:
    using LogicCriterion::LogicCriterion;

    unique_ptr<Criterion> unique_copy() override;
};

class LogicCriterionOR final : public LogicCriterion {
    [[nodiscard]] int check(int, int) const override;

public:
    using LogicCriterion::LogicCriterion;

    unique_ptr<Criterion> unique_copy() override;
};

#endif

//! @}
